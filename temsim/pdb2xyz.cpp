/*  pdb2xyz.cpp
    
    short program to convert PDB (protein data base) files
    to my xyz format file - NOT guaranteed to work on all PDB files

------------------------------------------------------------------------
Copyright 2009-2016 Earl J. Kirkland

This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

---------------------- NO WARRANTY ------------------
THIS PROGRAM IS PROVIDED AS-IS WITH ABSOLUTELY NO WARRANTY
OR GUARANTEE OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
INCLUDING BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
IN NO EVENT SHALL THE AUTHOR BE LIABLE
FOR DAMAGES RESULTING FROM THE USE OR INABILITY TO USE THIS
PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA
BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR
THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH
ANY OTHER PROGRAM).

------------------------------------------------------------------------
        
    started 2-nov-2009 E. Kirkland
    add carbon support with sub. from biocells.cpp 8-nov-2009
    convert to strings and streams 18-nov-2014 ejk
    convert to vector<> 24-nov-2014 ejk
    add rotation 6-may-2015 ejk
    add option to make just a-carbon film 4-sep-2016 ejk
    convert malloc2D() to vector<>  9-sep-2016 ejk
    change RMIN from 1.0 to 1.5 Ang (C-C bond) 10-sep-2016 ejk
    last updated 10-sep-2016 ejk
*/

#include <cstdlib>
#include <cmath>
#include <ctype.h>   // for toupper() 
#include <string.h>  // for strings
#include <time.h>

#include <string>
#include <iostream>  //  C++ stream IO
#include <sstream>   // string streams
#include <fstream>
#include <iomanip>   //  to format the output
#include <vector>

using namespace std;

#include "slicelib.hpp"

const int TRUE=1;
const int FALSE=0;

const double CWEIGHT=12.01115;  // molec. weight of carbon
const double CDENSITY=2.0;      // density in gm/cm^3 approx. for amorphous C
const double CFILL_FRACTION=0.9;// filling fraction for this density
const double NAV=6.0225e23;     // Avagadro's number (#/mole)
const double RMIN=1.4;          // min separation distance of C (in Angstroms)

const int OX=0, OY=1, OZ=2; // index for atomic coordinates
const int NVAL=3;

unsigned long iseed;        // random number generator seed

/*  subroutines define at the end of the file */
void newcoord( vector<double> &c, 
    const double xmax, const double ymax, const double zmax );
void insertcoord( vector<double> &ctest, vector<vector<double>> &coord, 
        const int ncoord, const int pos );
int testcoord( vector<double> &ctest, vector<vector<double>> &coord, 
    const int ncoord, int *pos, const double rmin );
int fillSolid( vector<vector<double>> &coord, const int ncoord, 
    double ax, double by, double cz, double rmin );


int main()
{
    string cline, filin, filout, symb;
    int i, ic, iz,  np, nlen, ntotal, ctop, lreadPDB;
    long ltime;
    double xpos, ypos, zpos, occ, mytime, x0, y0, z0, cr,sr,ct,st;
    double xmin, xmax, ymin, ymax, zmin, zmax, rotat, tilt, pi;
    double ax, by, cz, cthick, xoff, yoff, zoff, wobble, density;

    vector<int> znum;
    vector<double> xp, yp, zp, oc;

    ifstream fpin;
    ofstream fpout;
    
    /*  the following are the chemical symbols for the periodic table */
    char symbol[] = {
                " HHeLiBe B C N O FNeNaMgAlSi P SCl"
                "Ar KCaScTi VCrMnFeCoNiCuZnGaGeAsSeBr"
                "KrRbSr YZrNbMoTcRuRhPdAgCdInSnSbTe"
                " IXeCsBaLaCePrNdPmSmEuGdTbDyHoErTm"
                "YbLuHfTa WReOsIrPtAuHgTlPbBiPoAtRn"
                "FrRaAcThPa UNpPuAmCmBkCfEsFmMdNoLr"
        };

    /*------- what this is ---------- */
    
    cout << "pdb2xyz version dated 10-sep-2016" << endl;
    cout << "  convert PDB data file to xyz format\n" << endl;
    cout << "  (may not work on all PDB files)\n" << endl;

    //------- convert symbols to upper case for PDB comparison ---------- 
    
    nlen = (int) strlen( symbol );
    for( i=0; i<nlen; i++) symbol[i] = toupper( symbol[i] );

    //----------- get file names and open ---- 
    
    cout << "Type name of input file with PDB data ('none' for just a-carbon):" << endl;
    cin >> filin;
    if( filin.compare( "none" ) == 0  ) lreadPDB = FALSE;
        else lreadPDB = TRUE;
    if( TRUE == lreadPDB ) {
        fpin.open( filin );
        if( fpin.bad() ) { cout << "Can't open file "<< filin<<endl; exit(0); }
    }

    cout << "Type name of output file to get xyz data:" << endl;
    cin >>  filout;
    fpout.open( filout );
    if( fpout.bad() ) { cout << "Can't open file" << filout<<endl; exit(0); }

    if( TRUE == lreadPDB ) {
        cout << "Type azimuthal and polar rotation angles (in degrees):" << endl;
        cin >> rotat >> tilt;
        pi = 4.0 * atan( 1.0 );
        rotat = rotat * pi/180.0;
        tilt  = tilt * pi/180.0;
    }

    ///----------- get carbon support info ---- 
    if( TRUE == lreadPDB ) {
        cout << "Type thickness of carbon support (<0 to disable):" << endl;
        cin >> cthick;
        cthick = fabs( cthick);
        if( (cthick > 0.0) && (TRUE == lreadPDB) ) 
            ctop = askYN( "Do you want the carbon support on the entrance instead of exit");
    }

    if( FALSE == lreadPDB ) {
        cout << "Type size of amorphous carbon x,y,z (in Ang.):" << endl;
        cin >> ax >> by >> cthick;

        if((ax < 0) || (by < 0) || ( cthick < 0) ) {
            cout << "need > 0 if not reading PDB file, cannot continue...." << endl;
            return( EXIT_SUCCESS );
        }
    }

    //  get CPU time for fun
    mytime = cputim();

    if( TRUE == lreadPDB ) {
        //--------- read data from file in complicated format -------
        //   first get total range and number of atoms 
        np = 0; 
        do {
            getline( fpin, cline ); // read a whole line

            //  select lines begin. with ATOM or HETAM  with atom coord.
            if( ( cline.find( "ATOM") == 0 ) ||
                ( cline.find( "HETATM") == 0 )  )  {
            
                //----  read x,y,z coordinates
                istringstream isbuf( cline.substr(30) );
                isbuf >> xpos >> ypos >> zpos >> occ;

                if( 0 == np ) {
                    xmin = xmax = xpos;
                    ymin = ymax = ypos;
                    zmin = zmax = zpos;
                } else {
                    if( xpos < xmin ) xmin = xpos;  // coord. range
                    if( xpos > xmax ) xmax = xpos;
                    if( ypos < ymin ) ymin = ypos;
                    if( ypos > ymax ) ymax = ypos;
                    if( zpos < zmin ) zmin = zpos;
                    if( zpos > zmax ) zmax = zpos;
                }
                xp.push_back( xpos );   //  save coord.
                yp.push_back( ypos );
                zp.push_back( zpos );
                oc.push_back( occ );

                //  find atomic number
                symb= cline.substr(76,77);  // get chemical symbol
                for( i=0; i<nlen; i+=2) {
                    iz = 1 + i/2;
                    if( strncmp( &symbol[i], symb.c_str(), 2 ) == 0 ) break;
                }
                znum.push_back( iz );

                np++;
            }

        }  while( cline.find( "END") == string::npos );
        fpin.close( );
    
        cout <<  "Total number of atoms = " << np << endl;
        cout <<  "  with x range " << xmin << " to " << xmax << endl;
        cout <<  "   and y range " << ymin << " to " << ymax << endl;
        cout <<  "   and z range " << zmin << " to " << zmax << endl;

        if( ( fabs( rotat) > 1.0e-6 ) || ( fabs( tilt ) > 1.0e-6 ) ) {
            //-----  rotate to get better view perhaps (from slicview.cpp)
            //  Move to center of molecule and rotate
            xoff = 0.5F*( xmax + xmin );
            yoff = 0.5F*( ymax + ymin );
            zoff = 0.5F*( zmax + zmin );

            cout << "rotate about x,y,z= " << xpos << ", " << ypos << ", " << zpos << endl;

            /*  Calculate misc constants  */
            cr = cos( rotat );
            sr = sin( rotat );
            ct = cos( tilt );
            st = sin( tilt );

            for( i=0; i<np; i++) {
                xp[i] = xp[i] - xoff;   // translate to center
                yp[i] = yp[i] - yoff;
                zp[i] = zp[i] - zoff;

                x0 = xp[i];      //  Rotation
                y0 = yp[i];
                xp[i] =  xpos = cr*x0 - sr*y0;
                yp[i] =         sr*x0 + cr*y0;

                y0 = yp[i];      //  Tilt
                z0 = zp[i];
                yp[i] = ypos =  ct*y0 + st*z0;
                zp[i] = zpos = -st*y0 + ct*z0;
                if( (fabs( y0-yp[i]) > 0.1) || (fabs(z0-zp[i])>0.1 ) )

                //  find new range
                if( 0 == np ) {
                    xmin = xmax = xpos;
                    ymin = ymax = ypos;
                    zmin = zmax = zpos;
                } else {
                    if( xpos < xmin ) xmin = xpos;
                    if( xpos > xmax ) xmax = xpos;
                    if( ypos < ymin ) ymin = ypos;
                    if( ypos > ymax ) ymax = ypos;
                    if( zpos < zmin ) zmin = zpos;
                    if( zpos > zmax ) zmax = zpos;
                }
            }

            cout <<  "New range of atom coord. after rotation and tilt" << endl;
            cout <<  "   x range " << xmin << " to " << xmax << endl;
            cout <<  "   y range " << ymin << " to " << ymax << endl;
            cout <<  "   z range " << zmin << " to " << zmax << endl;

        }  // end rotate section

        //--------- write data to xyz file with offset -------

        //ax = 1.5*( xmax - xmin );
        //by = 1.5*( ymax - ymin );
        ax = 1.4*( xmax - xmin );
        by = 1.4*( ymax - ymin );
        if( by > ax ) ax = by;      // make it square to look right 
        else by = ax;
        if( cthick > 0.0 ) cz = (zmax - zmin) + cthick;
        else cz = zmax - zmin;
        fpout << "pdb2xyz translation of " << filin << endl;
        fpout << setw(16) << ax << setw(16) << by << setw(16) << cz << endl; // cell size 

        xoff = 0.5*ax - 0.5*(xmax+xmin);   // move molecule to the center 
        yoff = 0.5*ax - 0.5*(ymax+ymin);
        if( (cthick>0.0) && (ctop==1) ) zoff = -zmin + cthick;
        else zoff = -zmin;
    
        wobble = 0.0;   //  Debye-Waller factor 
        for( i=0; i<np; i++ ) {
            fpout << setw(5) << znum[i] << setw(14) << xp[i]+xoff << setw(14)
                << yp[i]+yoff << setw(14) << zp[i]+zoff << setw(14) << oc[i] 
                //<< zp[i]+yoff << setw(14) << yp[i]+zoff << setw(14) << oc[i] 
                << setw(14) << wobble << endl;
        }

    } else {  //  end if( TRUE == lreadPDB) ...

        //  write header for a-carbon only
        fpout << "pdb2xyz a-carbon" << endl;
        fpout << setw(16) << ax << setw(16) << by << setw(16) << cthick << endl; // cell size 
    }

    //----------- generate carbon support if requested -------------- 

    if( cthick > 0.0 ) {

        // ---- initialize random number generator seed ---- 
        cout << endl;
        ltime = (long) time( NULL );
        iseed = (unsigned) ltime;
        if( ltime == -1 ) {
            cout << "Type initial seed for random number generator:" << endl;
            cin >> iseed;
        } else {
            cout <<  "Random number seed initialized to " << iseed << endl;
        }
        cout << endl;

        
        cout << "generate random coord. for carbon support"  << endl;
        density = CFILL_FRACTION*NAV*CDENSITY*(1.0e-24)/CWEIGHT; //  # atoms/Angs^3 
        cout << "average density = "<< CFILL_FRACTION*CDENSITY << " gm/cm^3 = "
            << density << " atoms/Angstrom^3" << endl;
        cout << "minimum allowed separation = " << RMIN << " Angstroms" << endl;

        cout << "calculate coord. in a vol. of "<< ax << " x "<< by << 
            " x "<< cthick << " Angstroms" << endl;
        ntotal = (int) ( ax*by*cthick  * density );
        cout << "Total number of carbon atoms = " << ntotal << endl;

        vector<double> tmp( NVAL );
        vector< vector<double> > coord( ntotal, tmp );

        //  fill in random coord.
        ic = fillSolid( coord, ntotal, ax, by, cthick, RMIN );

        iz = 6;   //  atomic number of carbon 
        occ = 1.0;
        wobble = 0.0;
        if( ctop == 1 ) zoff = 0.0; //  on top  
        else zoff = zmax - zmin;    // on bottom 
        for( i=0; i<ic; i++) 
            fpout << setw(5) << iz << setw(14) << coord[i][0] << setw(14)
                << coord[i][1] << setw(14) << zoff + coord[i][2] << setw(14) << occ 
                << setw(14) << wobble << endl;
    }

    //-----------  write end of file mark -------------------------------- 
    iz = -1;
    fpout << setw(5) << iz << endl;  //  end of data 
    fpout.close( );

    // ------- echo CPU time just for fun ------------ 

    mytime = cputim() - mytime;
    cout << "\ntotal CPU time = " << mytime << " sec." << endl;

    return( EXIT_SUCCESS );

}  // end main 


/*------------------------ fillSolid() ------------------------*/
/*
    fill a vol. with an amorphous solid with random coord.
    in the range (0,0,0) to (ax,by,cz)

  coord[][] = to get list of existing sorted by z = coord[][OZ]
                dimensions ncoord x NVAL
  ncoord = number of coordinates to generate
  rmin = minimum separation distance

  assumed globals
     NVAL

*/
int fillSolid( vector< vector<double> > &coord, const int ncoord, 
    double ax, double by, double cz, double rmin )
{
    int ic, i;
    vector<double> ctest(NVAL);

    newcoord( coord[0], ax, by, cz );   // start coord. list
    ic = 1;
    do {
        // get new coordinate 
        newcoord( ctest, ax, by, cz );
        //cout << "new coord= ", << ctest[0] << ", " << 
        //           ctest[1] << ", " <<ctest[2] << endl;

        //  add it to the list if its OK - this also sorts
        if( testcoord( ctest, coord, ic, &i, RMIN ) == TRUE ) {
            insertcoord( ctest, coord, ic, i);
            ic++;
            if( ic%1000 == 0 ) cout << "\r ic= " << ic;
        }
    }  while( ic < ncoord );

    return ic;

}  // end fillSolid()

/*------------------------ insertcoord() ------------------------*/
/*
    insert a new coordinate in the list

  ctest[] = new coord to insert
  coord[][] = list of existing coordinates sorted by z = coord[][OZ]
  ncoord = number of coordinates
  pos = position to insert coord. at

  insert new coord at index pos and move all the rest down one

  assumed globals
     NVAL

     might use vector<> insert here sometime
*/
void insertcoord( vector<double> &ctest, vector<vector<double>> &coord, 
        const int ncoord, const int pos )
{
    long i, j;
    
    for( i=ncoord; i>pos; i--) {
        for( j=0; j<NVAL; j++)
            coord[i][j] = coord[i-1][j];
    }

    for( j=0; j<NVAL; j++)
            coord[pos][j] = ctest[j];

} // end insertcoord() 


/*------------------------ newcoord() ------------------------*/
/*
    generate a new random coordinate inside the required volume

    xmax, ymax, zmax = volume size
*/
void newcoord( vector<double> &c, 
    const double xmax, const double ymax, const double zmax )
{
    c[OX] = xmax * ranflat( &iseed );
    c[OY] = ymax * ranflat( &iseed );
    c[OZ] = zmax * ranflat( &iseed );
    return;

} // end newcoord() 

/*------------------------ testcoord() ------------------------*/
/*
    test the new coordinate to see if its too close to an existing
    coordinate (closer than RMIN)

    A straight search of the whole list is very slow (proportional to N^2)
    and was not pratical for a large set of coordinates (i.e. took
    an absurd amount of CPU time).

    This routines assumes that the list of existing coordinates is
    sorted wrt one coord (z in this case).  Once the position of the new
    coordinate is located in the list (using a binary search), then
    only the small range of coord. near this point need to be
    tested.  This makes the test dramatically faster.

  ctest[] = new coord to test
  coord[][] = list of existing coordinates sorted by z = coord[][OZ]
  ncoord = number of coordinates
  pos = returned position to insert coord.

  assumed globals
     OX, OY, OZ, RMIN

*/
int testcoord( vector<double> &ctest, vector< vector<double> > &coord, 
    const int ncoord, int *pos, const double rmin )
{
    /* this is the more -sophisticated version sorted by Z */
    long i, j, k;
    double d, dx, dy, dz, dz2, rmin2, z, range;

    /* ---- find postion of ctest[] in coord[][] using a binary search  --- 
         i should get the position to insert (between i and i+1)
         this assumes coord[][] is sorted by coord OZ
    */
//  printf("testcoord() top, ncoord= %d\n", ncoord );

    z = ctest[OZ];
    if( z <= coord[0][OZ] ) i = j = 0;
    else if( z >= coord[ncoord-1][OZ] ) i = j = ncoord;
    else { 
        i = 0;
        j = ncoord-1;
        do{ k = ( i + j ) / 2 ;
            if( z < coord[k][OZ] )  j = k;
            else if( z >=  coord[k][OZ] ) i = k;
        } while ( (j-i) > 1 );
    }

    /* now that we have the position of the new point
       we only have to explore within RMIN of this point */

    rmin2 = rmin*rmin;
    range = 2.0*rmin2;   /* add a little safety margin */
    k = j;
    while( k >= 0 ) {
            dx = ctest[OX] - coord[k][OX];
            dy = ctest[OY] - coord[k][OY];
            dz = ctest[OZ] - coord[k][OZ];
            dz2 = dz*dz;
            d = dx*dx + dy*dy + dz2;
            if( d <= rmin2 ) return( FALSE );
            if( dz2 > range ) break;
            k--;
    }
    k = j;
    while( k < ncoord ) {
            dx = ctest[OX] - coord[k][OX];
            dy = ctest[OY] - coord[k][OY];
            dz = ctest[OZ] - coord[k][OZ];
            dz2 = dz*dz;
            d = dx*dx + dy*dy + dz2;
            if( d <= rmin2 ) return( FALSE );
            if( dz2 > range ) break;
            k++;
    }

    /* if it gets to here then this is a good point */

    *pos = j;

    return( TRUE );

} // end testcoord() 

