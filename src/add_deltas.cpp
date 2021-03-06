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
#include "add_deltas.c"

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
    const size_t I = 1, O = 1;
    ifstream ifs1; ofstream ofs1;
    int8_t stdi1, stdo1, wo1;
    ioinfo i1, o1;
    uint32_t r, c;
    int N, dim;


    //Description
    string descr;
    descr += "Gets 1st order differences (deltas) of X along dim,\n";
    descr += "and appends them to X, so that Y has twice the size of X.\n";
    descr += "\n";
    descr += "The method is as used typically in ASR:\n";
    descr += "a linear ramp FIR from -N to N, with edges treated\n";
    descr += "as if using numpy.pad(A,N,'edge') (repeat edge sample).\n";
    descr += "\n";
    descr += "Use -n (--N) to give the number of samps before/after current samp.\n";
    descr += "That is, the FIR has length 2*N+1, with values -N to N.\n";
    descr += "The default is 2, as used in Kaldi, but 4 also appears typical.\n";
    descr += "\n";
    descr += "Output (Y) has the same data type and file format as X.\n";
    descr += "\n";
    descr += "Examples:\n";
    descr += "$ add_deltas X -o Y \n";
    descr += "$ add_deltas X > Y \n";
    descr += "$ cat X | add_deltas -n4 > Y \n";


    //Argtable
    int nerrs;
    struct arg_file  *a_fi = arg_filen(nullptr,nullptr,"<file>",I-1,I,"input file (X)");
    struct arg_int    *a_n = arg_intn("n","N","<uint>",0,1,"num samps such that FIR length is 2*N+1 [default=2]");
    struct arg_int    *a_d = arg_intn("d","dim","<uint>",0,1,"dimension along which to operate [default=0]");
    struct arg_file  *a_fo = arg_filen("o","ofile","<file>",0,O,"output file (Y)");
    struct arg_lit *a_help = arg_litn("h","help",0,1,"display this help and exit");
    struct arg_end  *a_end = arg_end(5);
    void *argtable[] = {a_fi, a_n, a_d, a_fo, a_help, a_end};
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
    wo1 = (stdo1 || a_fo->count>0);


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

    //Get N
    if (a_n->count==0) { N = 2; }
    else if (a_n->ival[0]<1) { cerr << progstr+": " << __LINE__ << errstr << "N must be positive" << endl; return 1; }
    else { N = a_n->ival[0]; }

    //Get dim
    if (a_d->count==0) { dim = 0; }
    else if (a_d->ival[0]<0) { cerr << progstr+": " << __LINE__ << errstr << "dim must be nonnegative" << endl; return 1; }
    else { dim = a_d->ival[0]; }
    if (dim!=0 && dim!=1) { cerr << progstr+": " << __LINE__ << errstr << "dim must be 0 or 1" << endl; return 1; }


    //Checks
    if (!i1.ismat()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) must be a matrix" << endl; return 1; }
    if (i1.isempty()) { cerr << progstr+": " << __LINE__ << errstr << "input 1 (X) found to be empty" << endl; return 1; }
    if (dim==0 && N>int(i1.R)) { cerr << progstr+": " << __LINE__ << errstr << "N must be <= nrows X for dim==0" << endl; return 1; }
    if (dim==1 && N>int(i1.C)) { cerr << progstr+": " << __LINE__ << errstr << "N must be <= ncols X for dim==1" << endl; return 1; }


    //Set output header info
    o1.F = i1.F; o1.T = i1.T;
    o1.R = (dim==0) ? i1.R : 2*i1.R;
    o1.C = (dim==1) ? i1.C : 2*i1.C;
    o1.S = i1.S; o1.H = i1.H;


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
    if (i1.T==1)
    {
        float *Y;
        try { Y = new float[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        if ((o1.iscolmajor() && dim==0) || (o1.isrowmajor() && dim==1))
        {
            try { ifs1.read(reinterpret_cast<char*>(Y),i1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        }
        else if (dim==0)
        {
            for (r=0u; r<i1.R; r++)
            {
                try { ifs1.read(reinterpret_cast<char*>(&Y[r*o1.C]),i1.sz()*i1.C); }
                catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X) at row " << r << endl; return 1; }
            }
        }
        else
        {
            for (c=0u; c<i1.C; c++)
            {
                try { ifs1.read(reinterpret_cast<char*>(&Y[c*o1.R]),i1.sz()*i1.R); }
                catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X) at col " << c << endl; return 1; }
            }
        }
        if (ov::add_deltas_s(Y,o1.iscolmajor(),int(o1.R),int(o1.C),dim,N)) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
        if (wo1)
        {
            try { ofs1.write(reinterpret_cast<char*>(Y),o1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem writing output file (Y)" << endl; return 1; }
        }
        delete[] Y;
    }
    else if (i1.T==2)
    {
        double *Y;
        try { Y = new double[o1.N()]; }
        catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem allocating for output file (Y)" << endl; return 1; }
        if ((o1.iscolmajor() && dim==0) || (o1.isrowmajor() && dim==1))
        {
            try { ifs1.read(reinterpret_cast<char*>(Y),i1.nbytes()); }
            catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X)" << endl; return 1; }
        }
        else if (dim==0)
        {
            for (r=0u; r<i1.R; r++)
            {
                try { ifs1.read(reinterpret_cast<char*>(&Y[r*o1.C]),i1.sz()*i1.C); }
                catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X) at row " << r << endl; return 1; }
            }
        }
        else
        {
            for (c=0u; c<i1.C; c++)
            {
                try { ifs1.read(reinterpret_cast<char*>(&Y[c*o1.R]),i1.sz()*i1.R); }
                catch (...) { cerr << progstr+": " << __LINE__ << errstr << "problem reading input file (X) at col " << c << endl; return 1; }
            }
        }
        if (ov::add_deltas_d(Y,o1.iscolmajor(),int(o1.R),int(o1.C),dim,N)) { cerr << progstr+": " << __LINE__ << errstr << "problem during function call" << endl; return 1; }
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

