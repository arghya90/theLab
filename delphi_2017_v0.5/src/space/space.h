#ifndef SPACE_H
#define SPACE_H

#include <complex>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <memory>
#include <valarray>
//#include <deque>

#include <cmath>      // std::abs

#include "../interface/interface_abstractmodule.h"
#include "../misc/misc_timer.h"
#include "../delphi/delphi_constants.h"
#include "../io/io.h"
#include "../interface/interface_datacontainer.h"
#include "space_templates.h"
#include "space_exceptions.h"
#include <stdio.h>
#include <stdlib.h>

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;
using namespace std;

class CDelphiSpace:virtual public IAbstractModule
{
private:                                            // In DATA CONTAINER
    shared_ptr<CTimer> pTimer;

    /*********************************************************************************************
    *                                                                                           *
    *              references to the variables obtained from the data container                 *
    *                                                                                           *
    ********************************************************************************************/

    //++++++++++++++ const references to read-only variables from data container +++++++++++++++//

    delphi_integer& iNatom;                           // natom
    const delphi_real&    fScale;                           // scale
    const delphi_integer& iGrid;                            // igrid (modified in setFocusBndy)
    const delphi_integer& iNObject;                         // nobject
    const delphi_real&    repsout;                           // repsout
    const delphi_real&    repsin;                           // repsin
	const delphi_real& fDebyeLength;
	const delphi_real& fEpsOut;
	const delphi_real& fEpsIn;

    //ARGO 14-FB,2016
    //To use for vdwSurf_file to be written as a CUBE file
    const SGrid<delphi_real>& fgBoxCenter;                   // oldmid
    //
    const bool& bUniformDiel;                       //
    //const bool& bVerbose;                           // removed
    const bool& bOnlyMol;                           //
    const bool& isolv;                           //
    const bool& irea;                           //
    const bool& logs;                           //
    const bool& lognl;                           //
    const bool& isen;                           //
    const bool& isch;                           //
    const bool& isite;                           //
    const bool& ibem;                           //
    const int& ibctyp;
    const bool& isitsf;
    const bool& bEpsOut;                        // only for the special case of CONVOLUTE

    //----- io file names ---  Required for CONVOLUTE EpsMap
    const string&  strEpsFile;
    //const bool& bDebug;                             //

    //const int& iTestGloble;                         //
    //const int& ibmx;                                //ibmx is removed

    const delphi_integer& iNMedia;
    const delphi_integer& numbmol;

    const delphi_integer& scrgfrm;
    //const delphi_integer& ndistr;
    const SGrid<delphi_real>& cOldMid;

    const delphi_real& fIonStrenth;
    const delphi_real& fExternRadius;
    const delphi_real& fRMax;
    const vector < delphi_real >& fRadPrb_v;
    const delphi_real* fRadPrb;

    //ARGO 
    // a new 2-element vector containing the minmax atomic indices from PDB
    // used for modifying Gaussian
    vector<delphi_integer>& minmax_AtomIndex_v;


    //const vector <CAtomPdb>& delphipdb;              //delphipdb

    const vector <string> & dataobject_v;




    //++++++++++++++++ reference to read-and-write variables from data container +++++++++++++++//

    delphi_integer ndistr;                         // ndistr

    delphi_integer& iBoundNum;                       // ibnum
    delphi_real& rdmx;                              // rdmx
    delphi_integer& nqass;
    delphi_integer& nqgrd;
    bool& iacs;                               //
    bool& isrf;
    SGrid<delphi_real>& cMin;
    SGrid<delphi_real>& cMax;
    delphi_real& qnet;
    delphi_real& qmin;
    delphi_real& qplus;
    SGrid<delphi_real>& cqmin;
    SGrid<delphi_real>& cqplus;
    SGrid<delphi_real>& acenter; //for focusing

    vector <delphi_real>& medeps;
    vector < SGrid<delphi_real> >& xn1_v;      //prgfgAtomCoordA
    vector < SGrid<delphi_real> >& xn2_v;      //prgfgAtomCoordG

    vector <bool>& bDebMap_v;             //idebmap


    vector < SGrid<delphi_integer> >& iEpsMap_v; //iepsmap
    vector < SGrid<delphi_real> >& fGepsMap_v; //gepsmap
    vector < SGrid<delphi_real> >& fGepsMap2_v; //gepsmap2
	vector < delphi_real >& fGDensityMap_v; //gaussian density map

    vector <delphi_integer>& iAtomMed_v;                 //iatmmed
    vector < SExtrema<delphi_real> >& sLimObject;       //limobject

    vector <SGrid <delphi_integer> >& ibgrd_v;               //ibgrd

    vector < SGrid <delphi_real> >& scspos_v;

    vector < SGrid <delphi_real> >& chgpos_v;

    vector <delphi_integer>& crgatn_v;
    vector <delphi_integer>& nqgrdtonqass_v;
    vector <delphi_real>& atmeps_v;
    vector< SGridValue<delphi_real> >& atmcrg_v;
    vector< SGridValue<delphi_real> >& chrgv2_v;

    vector < SGrid <delphi_real> >& scsnor_v;
    vector < delphi_integer >& atsurf_v;
    vector < delphi_integer >& atndx_v;

    vector <CAtomPdb>& delphipdb;              //delphipdb

    const float& cutoff;
    const float&  sigma;
    int&  inhomo;
    const float&  srfcut;
    const int&  iGaussian;

    //+++++++++++++++++++ NON-refereces: +++++++++++++++++++++++++++++++++++++++++++++
    struct delphipdb_struc
    {
        SGrid <delphi_real> xyz;
        delphi_real charge;
        delphi_real radius;
        string atom_info;  //ARGO
    };

    delphipdb_struc * sDelPhiPDB;




    //################# semi global variables in this class #########################

    //bool debug_space = true;
    bool debug_space ;
    delphi_integer iBoundNumsurf,extot,iall,lcb, mcb, ncb, lcb1, mcb1, ncb1;
    delphi_real radpmax,grdi, cbai, sideinter, sidemin;
    SGrid <delphi_real> mnxyz, xyzo, mxxyz;
    SGrid <delphi_integer> lmncb1, lmncb;
    SExtrema <delphi_integer> LimEps;
    delphi_integer extracrg;

	//delphi_real  fEpsOut;
	//delphi_real  fDebyeLength;
	//delphi_real  fEpsIn;

	delphi_real  test1,test2,test3;
	//delphi_real  test2;
	//delphi_real  test3;

    vector <delphi_integer>  iab1_v, iab2_v, icume, ast, cbn1_v, cbn2_v, cbal, icbn;
    vector <delphi_real> r0,r02,rs2;
    vector < SGrid <delphi_real> > expos;
    vector < SExtrema<delphi_real> > sLimGridUnit;


    vector < SGrid <delphi_real> > vert, vnorm, vnorm2;
    vector < SGrid <delphi_integer> > vindx;
    vector <delphi_integer> vtlen, vtlst, vtpnt;

    delphi_integer ** tmlst;

    //SGrid <delphi_real> * scsnor;
    //delphi_integer * atsurf;
    //delphi_integer * atndx;



    SGrid<delphi_real> * xn1;
    SGrid<delphi_real> * xn2;

    SGrid <delphi_integer> * ibgrd;               //ibgrd
    SGridValue <delphi_real> * atmcrg;
    SGridValue <delphi_real> * chrgv2;
    SGrid <delphi_real> * scspos;
    SGrid <delphi_real> * chgpos;

    SGrid <delphi_real> * scsnor;
    delphi_integer * atsurf;
    delphi_integer * atndx;



    delphi_integer * crgatn;
    delphi_integer * nqgrdtonqass;
    delphi_real * atmeps;


    SGrid <delphi_integer> *** egrid;
    //bool *** bDebMap;
    bool *** idebmap;


    //ARGO: defining a new variable to store bools for extended vdw surface for zeta
    bool *** zetaSurfMap;
    int& zetaOn;
    delphi_real& zetaDistance;
    vector <bool>& zetaSurfMap_v;			//zetaSurfMap
    const string&  strZetaPhiFile;			//zphinam
    vector < delphi_real > dMoment;
    vector < delphi_real > qMoment;
    vector <delphi_real>& surf_grid_coords_v;
    vector <delphi_real>& surf_grid_index_v;
    //

    //ARGO UA 2016, EPSILON SMOOTHING
    //EPS MODEL1
    int iConvolute;
    delphi_real *** ginit_rhomap;	//will store the initially assigned rho[in] values to obtain bubble surface
    delphi_real ***	cepsmap;		// will be assigned corresponding EPS values based on ginit_rhomap
    delphi_real ***  HRhomap;		// HS function map based on the binarization of cepsmap using that EPS cutoff. (Argo you know what EPS to set  the boundary).
    delphi_real& fksigma;			// the kernal sigma. Proabbly will be used only the FFT/iFFT module.
    vector<delphi_real>& fHRhomap_v;
    delphi_real& fhvsd_eps;			//the EPSILON value that will be used to draw the distinction between inside and out.
    //bool debug_convolute = false;
    bool debug_convolute ;
    //

    //Distance Matrix
    delphi_real ** distMatrix2;           //squared distances
    bool ** adjMatrix;               //adjecant matrix


    SGrid <delphi_integer> *** iepsmp;
    SGrid <delphi_real> *** gepsmp;
    SGrid <delphi_real> *** gepsmp2;

	/**
	* Gaussian Density Map array
	* Stores the density of atoms on each grid point
	* 
	*/
	delphi_real *** gDensityMapOnGridPoint;

    delphi_integer *** cbn1, *** cbn2, *** iab1, *** iab2;
    delphi_integer * iAtomMed;


    SGrid <delphi_real> sgrid_temp_real;
    SGrid <delphi_integer> sgrid_temp_int;
    SGrid <delphi_real> sgrid_rho_real;


/*
    template <class T> T *** getKey_Ptr3( vector <T> v_1d, const int& iRows, const int& iColumns, const int& iPages)
    {
        //if (!keyExists(strKey)) throw CInexistentKey(strKey);

        //vector<T> * nConstVectorPtr = any_cast< vector<T> >(&myData[strKey]);

        //if (nConstVectorPtr->size() != (unsigned int)iRows*iColumns*iPages) throw CUnmatchSize(strKey);
        cout << "in getkey_ptr3: " << endl;
        if (v_1d.size() != (unsigned int)iRows*iColumns*iPages) cerr << "error while getting ptr" << endl;

        T * nConstDataPtr = v_1d.data();

        T *** prg3D = new T ** [iRows];

        for (int i = 0; i < iRows; i++)
        {
            prg3D[i] = new T * [iColumns];

            for (int j = 0; j < iColumns; j++)
                prg3D[i][j] = &nConstDataPtr[i*iColumns*iPages+j*iPages];
        }

        return prg3D;
    }
*/

    //############### Functions in other files: ######################
    void epsmak();
    void setout();
    void setGaussian();
    void VdwToMs();
    void VdwToMs_piece(bool & , const delphi_integer& , const delphi_integer& , const delphi_integer& , const delphi_integer& , const delphi_real& , const SGrid <delphi_real>& , delphi_integer * , delphi_integer * , delphi_integer &  ,delphi_real & );
    void sas();
    void cube();
    void cubedata(delphi_real, delphi_real );
    void indverdata(delphi_real);
    void indver(delphi_integer);
    void sclbp();
    void msrf();
    void crgarr();


    //ARGO UA 2016 FFT
    void setConvolute();
    void Convolute ();
    CArray gaussianKernel( int );
    CArray fft( CArray& );
    CArray ifft( CArray & );

    //ARGO distance Matrix
    void CalculateDistMatrix();
    delphi_real calculateDistance2 ( SGrid <delphi_real> , SGrid <delphi_real> );


    /**
       * function to remove spaces from a string
       * @param strLine input string
       * @return converted string without spaces
       */
    string removeSpace(const string& strLine);




    SGrid <int> int_coord( const int& a, const int& b,  const int& c);

    SGrid <float> coord( const float& a,  const float& b,  const float& c);

    SGrid <int> Float2Int( const SGrid <float>& a );

    int Float2Int(const float& a);

    SGrid <float> Int2Float( const SGrid <int>& a );

    float Int2Float(const int& a);

    shared_ptr<IDataContainer> test_pdc;



public:
    SGrid <delphi_integer> *** iEpsMap;
    CDelphiSpace(shared_ptr<IDataContainer> pdc,shared_ptr<CTimer> pt):
/*********************************************************************************************
*                                                                                           *
*              references to the variables obtained from the data container                 *
*                                                                                           *
********************************************************************************************/

        //++++++++++++++ const references to read-only variables from data container +++++++++++++++//
        IAbstractModule(pdc),
        pTimer(pt),

        iNatom (pdc->getKey_Ref<delphi_integer>("natom")),
        fScale (pdc->getKey_constRef<delphi_real>("scale")),
        iGrid (pdc->getKey_constRef<delphi_integer>("igrid")),
        iNObject (pdc->getKey_constRef<delphi_integer>("nobject")),
        repsout (pdc->getKey_constRef<delphi_real>("repsout")),
        repsin (pdc->getKey_constRef<delphi_real>("repsin")),
        //ndistr (pdc->getKey_constRef<delphi_integer>("ndistr")),
        bUniformDiel (pdc->getKey_constRef<bool>("uniformdiel")),
        bOnlyMol (pdc->getKey_constRef<bool>("ionlymol")),
        isolv (pdc->getKey_constRef<bool>("isolv")),
        irea (pdc->getKey_constRef<bool>("irea")),
        logs (pdc->getKey_constRef<bool>("logs")),
        lognl (pdc->getKey_constRef<bool>("lognl")),
        isen (pdc->getKey_constRef<bool>("isen")),
		fDebyeLength(pdc->getKey_constRef<delphi_real>("deblen")),
        isch (pdc->getKey_constRef<bool>("isch")),
        fEpsOut(pdc->getKey_constRef<delphi_real>("epsout")),
		fEpsIn(pdc->getKey_constRef<delphi_real>("epsin")),

        //Argo: for EPSMAP in conjunction with CONVOLUTION
        bEpsOut(pdc->getKey_constRef<bool>("epswrt")),
        strEpsFile(pdc->getKey_constRef<string>("epsnam")),
        //ARGO-Putting in the reference for fgBoxCenter
		fgBoxCenter(pdc->getKey_constRef< SGrid<delphi_real> >("oldmid")),
		    //
        isite (pdc->getKey_constRef<bool>("isite")),
        ibem (pdc->getKey_constRef<bool>("ibem")),
        ibctyp (pdc->getKey_constRef<int>("ibctyp")),
        isitsf (pdc->getKey_constRef<bool>("isitsf")),

        // Lin Li : Gaussian:
        cutoff (pdc->getKey_constRef<float>("cutoff")),
        sigma (pdc->getKey_constRef<float>("sigma")),
        inhomo (pdc->getKey_Ref<int>("inhomo")),		//to be used for CONVOLUTE also
        srfcut (pdc->getKey_constRef<float>("srfcut")),
        iGaussian (pdc->getKey_constRef<int>("gaussian")),
        minmax_AtomIndex_v (pdc->getKey_Ref< vector<delphi_integer> >("MinMaxAtomIdx")), //ARGO: modified Gaussian


        //bDebug (pdc->getKey_constRef<bool>("debug")),

        //iTestGloble (pdc->getKey_constRef<int>(" ")),
        iNMedia (pdc->getKey_constRef<delphi_integer>("nmedia")),
        numbmol (pdc->getKey_constRef<delphi_integer>("numbmol")),
        scrgfrm (pdc->getKey_constRef<delphi_integer>("scrgfrm")),

        //ndistr (pdc->getKey_constRef<delphi_integer>("ndistr")),

        //fRMid (pdc->getKey_constRef<delphi_real>("rmid")),
        cOldMid(pdc->getKey_constRef< SGrid<delphi_real> >("oldmid")),
        fIonStrenth (pdc->getKey_constRef<delphi_real>("rionst")),
        fExternRadius (pdc->getKey_constRef<delphi_real>("exrad")),
        fRMax (pdc->getKey_constRef<delphi_real>("rdmx")),
        fRadPrb_v (pdc->getKey_constRef<vector <delphi_real> >("radprb")),

        dataobject_v(pdc->getKey_constRef< vector<string> >("dataobject")),

        //++++++++++++++++ reference to read-and-write variables from data container +++++++++++++++//

        delphipdb (pdc->getKey_Ref<vector <CAtomPdb> >("delphipdb")),
        iBoundNum (pdc->getKey_Ref<delphi_integer>("ibnum")),
        rdmx (pdc->getKey_Ref<delphi_real>("rdmx")),
        nqass (pdc->getKey_Ref<delphi_integer>("nqass")),
        nqgrd (pdc->getKey_Ref<delphi_integer>("nqgrd")),
        iacs(pdc->getKey_Ref< bool >("iacs")),
        isrf(pdc->getKey_Ref< bool >("isrf")),
        cMin(pdc->getKey_Ref< SGrid<delphi_real> >("cmin")),
        cMax(pdc->getKey_Ref< SGrid<delphi_real> >("cmax")),
        qnet(pdc->getKey_Ref< delphi_real >("qnet")),
        qmin(pdc->getKey_Ref< delphi_real >("qmin")),
        qplus(pdc->getKey_Ref< delphi_real >("qplus")),
        cqmin(pdc->getKey_Ref< SGrid<delphi_real> >("cqmin")),
        cqplus(pdc->getKey_Ref< SGrid<delphi_real> >("cqplus")),
        acenter(pdc->getKey_Ref< SGrid<delphi_real> >("acent")),
        medeps(pdc->getKey_Ref < vector < delphi_real > >("medeps")),
        xn1_v(pdc->getKey_Ref< vector< SGrid<delphi_real> > >("xn1")),
        xn2_v(pdc->getKey_Ref< vector< SGrid<delphi_real> > >("xn2")),
        bDebMap_v( pdc->getKey_Ref< vector< bool > >("idebmap")),

        //ARGO: doing the 'pdc' thing just like its done with idebmap
        zetaSurfMap_v( pdc->getKey_Ref< vector< bool > >("zetaSurfMap")),
        zetaOn( pdc->getKey_Ref<int>("zetaOn")),
        zetaDistance( pdc->getKey_Ref<delphi_real>("zetaDistance")),
        strZetaPhiFile(pdc->getKey_constRef<string>("zphinam")),
        dMoment( pdc->getKey_Ref< vector<delphi_real> >("dMoment")),
        qMoment( pdc->getKey_Ref< vector<delphi_real> >("qMoment")),
        surf_grid_coords_v( pdc->getKey_Ref< vector<delphi_real> >("surf_grid_coords") ),
        surf_grid_index_v( pdc->getKey_Ref< vector<delphi_real> >("surf_grid_index") ),
    
        //

        //ARGO: doing the 'pdc' thing just like others
        iConvolute(pdc->getKey_Ref<int> ("convolute")),
        fksigma( pdc->getKey_Ref<delphi_real>("kernel_sigma")),
        fHRhomap_v( pdc->getKey_Ref< vector<delphi_real> >("rho_map")),
        fhvsd_eps(pdc->getKey_Ref<delphi_real>("hvsd_epsilon")),
        //

        iEpsMap_v( pdc->getKey_Ref< vector< SGrid<delphi_integer> > > ("iepsmp")),
        fGepsMap_v( pdc->getKey_Ref< vector< SGrid<delphi_real> > > ("gepsmp")),
        fGepsMap2_v( pdc->getKey_Ref< vector< SGrid<delphi_real> > > ("gepsmp2")), 
		fGDensityMap_v(pdc->getKey_Ref< vector< delphi_real> >("gdensity")),
        iAtomMed_v(pdc->getKey_Ref< vector<delphi_integer> >("iatmmed")),
        sLimObject(pdc->getKey_Ref< vector < SExtrema<delphi_real> > >("limobject")),
        ibgrd_v(pdc->getKey_Ref< vector< SGrid<delphi_integer> > >("ibgrd")),
        scspos_v(pdc->getKey_Ref< vector< SGrid<delphi_real> > >("scspos")),
        chgpos_v(pdc->getKey_Ref< vector< SGrid<delphi_real> > >("chgpos")),
        crgatn_v(pdc->getKey_Ref< vector< delphi_integer > >("crgatn")),
        nqgrdtonqass_v(pdc->getKey_Ref< vector< delphi_integer > >("nqgrdtonqass")),
        atmeps_v(pdc->getKey_Ref< vector< delphi_real > >("atmeps")),
        atmcrg_v(pdc->getKey_Ref< vector< SGridValue<delphi_real> > >("atmcrg")),
        chrgv2_v(pdc->getKey_Ref< vector< SGridValue<delphi_real> > >("chrgv2")),

        scsnor_v(pdc->getKey_Ref< vector< SGrid<delphi_real> > >("scsnor")),
        atsurf_v(pdc->getKey_Ref< vector< delphi_integer > >("atsurf")),
        atndx_v(pdc->getKey_Ref< vector< delphi_integer > >("atndx"))

    {
        //bDebMap_v=pdc->getKey_Ref< vector< bool > >("idebmap");
        bDebMap_v.assign(iGrid*iGrid*iGrid, true);

        //ARGO
        zetaSurfMap_v.assign(iGrid*iGrid*iGrid, true);


        //iEpsMap_v.assign(iGrid*iGrid*iGrid, {0,0,0});

        iAtomMed=&iAtomMed_v[-1];


        sDelPhiPDB = new delphipdb_struc [iNatom+1];
        //bDebMap=get_pt3d <bool> (iGrid,iGrid,iGrid);
        //get_pt3d <bool> (bDebMap,iGrid,iGrid,iGrid);

        //egrid = get_pt3d <SGrid <delphi_integer> > (iGrid,iGrid,iGrid);


        ndistr = 0;


        test_pdc=pdc;

        /*
        delphi_integer ** tmlst;
        SGrid <delphi_integer> *** egrid;
        bool *** idebmap;
        SGrid <delphi_integer> *** iepsmp;
        delphi_integer *** cbn1, *** cbn2, *** iab1, *** iab2;
        */

		//fEpsOuts = pdc->getKey_Val<delphi_real>("epsout");
		//fDebyeLengths = pdc->getKey_Val<delphi_real>("deblen");
		//fEpsIns = pdc->getKey_Val<delphi_real>("epsin");


		
        // initialize all the pointers to be NULL:
        tmlst=NULL;
        egrid=NULL;
        idebmap=NULL;

        //ARGO
        zetaSurfMap=NULL;
        ginit_rhomap=NULL;
        cepsmap=NULL;
        HRhomap=NULL;
        //

        //ARGO
        distMatrix2 = NULL;
        adjMatrix = NULL;

        iepsmp=NULL;
        gepsmp=NULL;
        gepsmp2=NULL;
		gDensityMapOnGridPoint=NULL;

        iab1=NULL;
        iab2=NULL;
        cbn1=NULL;
        cbn2=NULL;

    };


    ~CDelphiSpace() {
        delete [] sDelPhiPDB ;
        //free_pt3d <bool> (idebmap,iGrid+1,iGrid+1,iGrid+1);
        //free_pt3d <SGrid <delphi_integer> > (iEpsMap,iGrid+1,iGrid+1,iGrid+1);

    };


    virtual void validateInput();

    virtual void run();

};

#endif // SPACE_H
