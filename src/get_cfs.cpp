//@author Erik Edwards
//@date 2019-2020


#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string>
#include <cstring>
#include <valarray>
#include <complex>
#include <unordered_map>
#include <argtable2.h>
#include "/home/erik/codee/cmli/cmli.hpp"
#include "/home/erik/codee/openvoice/openvoice.h"


int main(int argc, char *argv[])
{
    using namespace std;


    //Declarations
    int ret = 0;
    const string errstr = ": \033[1;31merror:\033[0m ";
    const string warstr = ": \033[1;35mwarning:\033[0m ";
    const string progstr(__FILE__,string(__FILE__).find_last_of("/")+1,strlen(__FILE__)-string(__FILE__).find_last_of("/")-5);
    const valarray<uint8_t> oktypes = {1,2};
    const size_t O = 1;
    ofstream ofs1;
    int8_t stdo1, wo1;
    ioinfo o1;
    int dim, B;
    char edges;
    string fscale;
    string::size_type c;
    double lofreq, hifreq;


    //Description
    string descr;
    descr += "Gets B center frequencies (B is number of frequency bands),\n";
    descr += "ranging from lofreq to hifreq, in units of Hz, but\n";
    descr += "equally spaced on the specified frequency scale.\n";
    descr += "\n";
    descr += "The available frequency scales are: \n";
    descr += "'bark', 'cochlea', 'erb', 'hz', 'mel', 'midi', 'octave', 'piano', \n";
    descr += "'ihcs' (inner hair cells), and 'sgcs' (spiral ganglion cells).\n";
    descr += "\n";
    descr += "Use -b (--ncfs) to give the number of cfs (B).\n";
    descr += "\n";
    descr += "Use -d (--dim) to give the nonsingleton dim of the output vec.\n";
    descr += "If d=0, then Y is a column vector [default].\n";
    descr += "If d=1, then Y is a row vector.\n";
    descr += "\n";
    descr += "Include -e (--edges) to interpret lofreq and hifreq as the.\n";
    descr += "lower and upper edges of the first and last bands, respectively.\n";
    descr += "The default is to interpret these as frequency-band centers,\n";
    descr += "such that cfs[0] = lofreq and cfs[B-1] = hifreq.\n";
    descr += "\n";
    descr += "Since this is a generating function (no inputs),\n";
    descr += "the output data type and file format can be specified by\n";
    descr += "-t and -f, respectively (these are the usual CMLI opts).\n";
    descr += "\n";
    descr += "Examples:\n";
    descr += "$ get_cfs -m'mel' -b40 -o Y \n";
    descr += "$ get_cfs -m'mel' -u10000 -b40 -e > Y \n";
    descr += "$ get_cfs -m'bark' -l20 -u20000 -b64 -d1 -t1 -f101 > Y \n";


    //Argtable
    int nerrs;
    struct arg_int    *a_b = arg_intn("b","ncfs","<uint>",0,1,"number of frequency bands (cfs) [default=23]");
    struct arg_str   *a_sc = arg_strn("m","fscale","<str>",0,1,"frequency scale [default='mel']");
    struct arg_dbl  *a_lof = arg_dbln("l","lofreq","<uint>",0,1,"lowest frequency (cfs[0]) in Hz [default=0]");
    struct arg_dbl  *a_hif = arg_dbln("u","hifreq","<uint>",0,1,"upper (highest) frequency (cfs[B-1]) in Hz [default=5000]");
    struct arg_lit  *a_edg = arg_litn("e","edges",0,1,"include to use lofreq and hifreq as band edges");
    struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"nonsingleton dimension [default=0 -> col vec]");
    struct arg_int *a_otyp = arg_intn("t","type","<uint>",0,1,"output data type [default=2 -> double]");
    struct arg_int *a_ofmt = arg_intn("f","fmt","<uint>",0,1,"output file format [default=102 -> colmajor]");
    struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");
    struct arg_lit *a_help = arg_litn("h","help",0,1,"display this help and exit");
    struct arg_end  *a_end = arg_end(5);
    void *argtable[] = {a_b, a_sc, a_lof, a_hif, a_edg, a_d, a_otyp, a_ofmt, a_fo, a_help, a_end};
    if (arg_nullcheck(argtable)!=0) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating argtable" << endl; return 1; }
    nerrs = arg_parse(argc, argv, argtable);
    if (a_help->count>0)
    {
        cout << "Usage: " << progstr; arg_print_syntax(stdout, argtable, "\n");
        cout << endl; arg_print_glossary(stdout, argtable, "  %-25s %s\n");
        cout << endl << descr; return 1;
    }
    if (nerrs>0) { arg_print_errors(stderr,a_end,(progstr+": "+to_string(__LINE__)+errstr).c_str()); return 1; }


    //Check stdout
    if (a_fo->count>0) { stdo1 = (strlen(a_fo->filename[0])==0 || strcmp(a_fo->filename[0],"-")==0); }
    else { stdo1 = (!isatty(fileno(stdout))); }
    wo1 = (stdo1 || a_fo->count>0);


    //Get options

    //Get fscale
    if (a_sc->count==0) { fscale = "mel"; }
    else
    {
    	try { fscale = string(a_sc->sval[0]); }
    	catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem getting string for frequency scale" << endl; return 1; }
    }
    for (c=0; c<fscale.size(); c++) { fscale[c] = char(tolower(fscale[c])); }

    //Get lofreq
    if (a_lof->count==0) { lofreq = 0.0; }
    else if (a_lof->dval[0]<0.0) { cerr << progstr+": " << __LINE__ << errstr << "lofreq must be nonnegative" << endl; return 1; }
    else { lofreq = a_lof->dval[0]; }

    //Get hifreq
    if (a_hif->count==0) { hifreq = 1.0; }
    else if (a_hif->dval[0]<=0.0) { cerr << progstr+": " << __LINE__ << errstr << "hifreq must be positive" << endl; return 1; }
    else { hifreq = a_hif->dval[0]; }
    if (hifreq<=lofreq) { cerr << progstr+": " << __LINE__ << errstr << "hifreq must be > lofreq" << endl; return 1; }
    if (hifreq>=100000.0) { cerr << progstr+": " << __LINE__ << errstr << "hifreq must be < 100 kHz" << endl; return 1; }

    //Get B
    if (a_b->count==0) { B = 23; }
    else if (a_b->ival[0]<2) { cerr << progstr+": " << __LINE__ << errstr << "B must be a positive int > 1" << endl; return 1; }
    else { B = a_b->ival[0]; }

    //Get edges
    edges = (a_edg->count>0);

    //Get dim
    if (a_d->count==0) { dim = 0; }
    else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
    else { dim = a_d->ival[0]; }
    if (dim>3) { cerr << progstr+": " << __LINE__ << errstr << "dim must be in {0,1,2,3}" << endl; return 1; }

    //Get o1.F
    if (a_ofmt->count==0) { o1.F = 102; }
    else if (a_ofmt->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be nonnegative" << endl; return 1; }
    else if (a_ofmt->ival[0]>255) { cerr << progstr+": " << __LINE__ << errstr << "output file format must be < 256" << endl; return 1; }
    else { o1.F = uint8_t(a_ofmt->ival[0]); }

    //Get o1.T
    if (a_otyp->count==0) { o1.T = 2; }
    else if (a_otyp->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "data type must be positive int" << endl; return 1; }
    else { o1.T = uint8_t(a_otyp->ival[0]); }
    if ((o1.T==oktypes).sum()==0)
    {
        cerr << progstr+": " << __LINE__ << errstr << "output data type must be in " << "{";
        for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1]) ? "}" : ","); }
        cerr << endl; return 1;
    }


    //Set output header info
    o1.R = (dim==0) ? uint32_t(B) : 1u;
    o1.C = (dim==1) ? uint32_t(B) : 1u;
    o1.S = (dim==2) ? uint32_t(B) : 1u;
    o1.H = (dim==3) ? uint32_t(B) : 1u;


    //Open output
    if (wo1)
    {
        if (stdo1) { ofs1.copyfmt(cout); ofs1.basic_ios<char>::rdbuf(cout.rdbuf()); } else { ofs1.open(a_fo->filename[0]); }
        if (!ofs1) { cerr << progstr+": " << __LINE__ << errstr << "problem opening output file 1" << endl; return 1; }
    }


    //Write output header
    if (wo1 && !write_output_header(ofs1,o1)) { cerr << progstr+": " << __LINE__ << errstr << "problem writing header for output file 1" << endl; return 1; }


    //Other prep


    //Process
    if (o1.T==1)
    {
        float *Y;
        try { Y = new float[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        if (ov::get_cfs_s(Y,B,float(lofreq),float(hifreq),fscale.c_str(),edges)) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] Y;
    }
    else if (o1.T==2)
    {
        double *Y;
        try { Y = new double[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        if (ov::get_cfs_d(Y,B,double(lofreq),double(hifreq),fscale.c_str(),edges)) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] Y;
    }
    else
    {
        cerr << progstr+": " << __LINE__ << errstr << "data type not supported" << endl; return 1;
    }
    

    //Exit
    return ret;
}

