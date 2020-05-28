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
#include "sig2ar_burg.c"

#ifdef I
#undef I
#endif


int main(int argc, char *argv[])
{
    using namespace std;


    //Declarations
    int ret = 0;
    const string errstr = ": \033[1;31merror:\033[0m ";
    const string warstr = ": \033[1;35mwarning:\033[0m ";
    const string progstr(__FILE__,string(__FILE__).find_last_of("/")+1,strlen(__FILE__)-string(__FILE__).find_last_of("/")-5);
    const valarray<uint8_t> oktypes = {1,2};
    const size_t I = 1, O = 2;
    ifstream ifs1; ofstream ofs1, ofs2;
    int8_t stdi1, stdo1, stdo2, wo1, wo2;
    ioinfo i1, o1, o2;
    int dim, P;
    char m0;


    //Description
    string descr;
    descr += "Does linear prediction for each row or col of X,\n";
    descr += "using the Burg autoregressive (AR) method,\n";
    descr += "which uses time-domain forward and backward predictions.\n";
    descr += "\n";
    descr += "Use -p (--P) to specify the number of LP coefficients.\n";
    descr += "\n";
    descr += "Use -d (--dim) to give the dimension along which to operate.\n";
    descr += "Default is 0 (along cols), unless X is a row vector.\n";
    descr += "\n";
    descr += "Include -m (--mean0) to zero the mean of each row or col of X [default=no].\n";
    descr += "This is done before computing AR coeffs (works along same dim as above).\n";
    descr += "\n";
    descr += "If dim==0, then Y has size P x C.\n";
    descr += "If dim==1, then Y has size R x P.\n";
    descr += "\n";
    descr += "Examples:\n";
    descr += "$ sig2ar_burg -p3 X -o Y -o V \n";
    descr += "$ sig2ar_burg -d1 -p5 X -o Y -o V \n";
    descr += "$ cat X | sig2ar_burg -p7 -m > Y \n";


    //Argtable
    int nerrs;
    struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
    struct arg_int    *a_p = arg_intn("p","P","<uint>",0,1,"number of LP coeffs [default=1]");
    struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to operate [default=0]");
    struct arg_lit   *a_m0 = arg_litn("m","mean0",0,1,"zero mean before computing [default=no]");
    struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output files (Y,V)");
    struct arg_lit *a_help = arg_litn("h","help",0,1,"display this help and exit");
    struct arg_end  *a_end = arg_end(5);
    void *argtable[] = {a_fi, a_p, a_d, a_m0, a_fo, a_help, a_end};
    if (arg_nullcheck(argtable)!=0) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating argtable" << endl; return 1; }
    nerrs = arg_parse(argc, argv, argtable);
    if (a_help->count>0)
    {
        cout << "Usage: " << progstr; arg_print_syntax(stdout, argtable, "\n");
        cout << endl; arg_print_glossary(stdout, argtable, "  %-25s %s\n");
        cout << endl << descr; return 1;
    }
    if (nerrs>0) { arg_print_errors(stderr,a_end,(progstr+": "+to_string(__LINE__)+errstr).c_str()); return 1; }


    //Check stdin
    stdi1 = (a_fi->count==0 || strlen(a_fi->filename[0])==0 || strcmp(a_fi->filename[0],"-")==0);
    if (stdi1>0 && isatty(fileno(stdin))) { cerr << progstr+": " << __LINE__ << errstr << "no stdin detected" << endl; return 1; }


    //Check stdout
    if (a_fo->count>0) { stdo1 = (strlen(a_fo->filename[0])==0 || strcmp(a_fo->filename[0],"-")==0); }
    else { stdo1 = (!isatty(fileno(stdout))); }
    if (a_fo->count>1) { stdo2 = (strlen(a_fo->filename[1])==0 || strcmp(a_fo->filename[1],"-")==0); }
    else { stdo2 = (!isatty(fileno(stdout)) && a_fo->count==1 && stdo1==0); }
    if (stdo1+stdo2>1) { cerr << progstr+": " << __LINE__ << errstr << "can only use stdout for one output" << endl; return 1; }
    wo1 = (stdo1 || a_fo->count>0); wo2 = (stdo2 || a_fo->count>1);


    //Open input
    if (stdi1) { ifs1.copyfmt(cin); ifs1.basic_ios<char>::rdbuf(cin.rdbuf()); } else { ifs1.open(a_fi->filename[0]); }
    if (!ifs1) { cerr << progstr+": " << __LINE__ << errstr << "problem opening input file" << endl; return 1; }


    //Read input header
    if (!read_input_header(ifs1,i1)) { cerr << progstr+": " << __LINE__ << errstr << "problem reading header for input file" << endl; return 1; }
    if ((i1.T==oktypes).sum()==0)
    {
        cerr << progstr+": " << __LINE__ << errstr << "input data type must be in " << "{";
        for (auto o : oktypes) { cerr << int(o) << ((o==oktypes[oktypes.size()-1]) ? "}" : ","); }
        cerr << endl; return 1;
    }


    //Get options

    //Get dim
    if (a_d->count==0) { dim = (i1.R==1u) ? 1 : 0; }
    else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
    else { dim = a_d->ival[0]; }
    if (dim!=0 && dim!=1) { cerr << progstr+": " << __LINE__ << errstr << "dim must be 0 or 1" << endl; return 1; }

    //Get L
    if (a_p->count==0) { P = 1; }
    else if (a_p->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "P must be positive" << endl; return 1; }
    else { P = a_p->ival[0]; }

    //Get m0
    m0 = (a_m0->count>0);


    //Checks
    if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input must be 1D or 2D" << endl; return 1; }
    if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input (X) found to be empty" << endl; return 1; }


    //Set output header infos
    o1.F = o2.F = i1.F;
    o1.T = o2.T = i1.T;
    o1.R = (dim==0) ? uint32_t(P) : i1.R;
    o1.C = (dim==1) ? uint32_t(P) : i1.C;
    o2.R = (dim==0) ? 1u : i1.R;
    o2.C = (dim==1) ? 1u : i1.C;
    o1.S = o2.S = i1.S;
    o1.H = o2.H = i1.H;


    //Open outputs
    if (wo1)
    {
        if (stdo1) { ofs1.copyfmt(cout); ofs1.basic_ios<char>::rdbuf(cout.rdbuf()); } else { ofs1.open(a_fo->filename[0]); }
        if (!ofs1) { cerr << progstr+": " << __LINE__ << errstr << "problem opening output file 1" << endl; return 1; }
    }
    if (wo2)
    {
        if (stdo2) { ofs2.copyfmt(cout); ofs2.basic_ios<char>::rdbuf(cout.rdbuf()); } else { ofs2.open(a_fo->filename[1]); }
        if (!ofs2) { cerr << progstr+": " << __LINE__ << errstr << "problem opening output file 2" << endl; return 1; }
    }


    //Write output headers
    if (wo1 && !write_output_header(ofs1,o1)) { cerr << progstr+": " << __LINE__ << errstr << "problem writing header for output file 1" << endl; return 1; }
    if (wo2 && !write_output_header(ofs2,o2)) { cerr << progstr+": " << __LINE__ << errstr << "problem writing header for output file 2" << endl; return 1; }


    //Other prep


    //Process
    if (i1.T==1)
    {
        float *X, *Y, *V;
        try { X = new float[i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
        try { Y = new float[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 1 (Y)" << endl; return 1; }
        try { V = new float[o2.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 2 (V)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        if (ov::sig2ar_burg_s(Y,V,X,i1.iscolmajor(),int(i1.R),int(i1.C),dim,P,m0)) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        if (wo2)
        {
            try { ofs2.write(reinterpret_cast<char*>(V),o2.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (V)" << endl; return 1; }
        }
        delete[] X; delete[] Y; delete[] V;
    }
    else if (i1.T==2)
    {
        double *X, *Y, *V;
        try { X = new double[i1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for input file (X)" << endl; return 1; }
        try { Y = new double[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 1 (Y)" << endl; return 1; }
        try { V = new double[o2.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file 2 (V)" << endl; return 1; }
        try { ifs1.read(reinterpret_cast<char*>(X),i1.nbytes()); }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        if (ov::sig2ar_burg_d(Y,V,X,i1.iscolmajor(),int(i1.R),int(i1.C),dim,P,m0)) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        if (wo2)
        {
            try { ofs2.write(reinterpret_cast<char*>(V),o2.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (V)" << endl; return 1; }
        }
        delete[] X; delete[] Y; delete[] V;
    }
    else
    {
        cerr << progstr+": " << __LINE__ << errstr << "data type not supported" << endl; return 1;
    }
    

    //Exit
    return ret;
}

