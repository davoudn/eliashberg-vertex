#include<complex>
#include<vector>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include "matrix.h"
#include <time.h>
#include <fftw3.h>
#include <tuple>
//#define  KSUM
const double kelvin2mev=0.0;
const double ev2mev=0.0;
const std::complex<double> cone(1.0,0.0);
const std::complex<double> czero(0.0,0.0);
const std::complex<double> imath(0.0,1.0);
const double zero(0.0);
const double one(1.0);
const double pi = acos(-1.0);
const double invsqrtpi =1.0/1.7724538509055160; // std::pow(pi,-0.5);
const std::complex<double> cpi(pi,0.0);
const std::complex<float> fczero(0.0,0.0);


const double bohr     = 0.52917721092 ; 
const double ryd2mev  = 13605.6981  ;
const double ryd2ev   = 13.6056981  ;
const double rydcm1   = 13.6056981   * 8065.541 ;
const double bohr2ang = 0.52917721092  ;
const double ev2cmm1  = 8065.541  ;
const double kelvin2eV= 8.6173427909e-05;

//
template <typename T> int sgn(T val, T comp) {
    return (comp < val) - (val < -comp);
}

class Delta {
public:
Delta(){}
~Delta(){}
public:
double operator() (double x) { return invsqrtpi*exp(-pow(x,2.0));}

};

struct BOSON{};
struct FERMI{};
template <typename X>
struct Mesh;

template <>
struct Mesh<FERMI> {
using type = double;
Mesh(){}
Mesh(double T_,int size_): T(T_),mesh_val(size_),size__(size_){
    for (int i = -size_/2; i < size_/2; ++i ) mesh_val[i+size_/2] = pi*(2.0*i+1)*T;
//      std::cout << size_ << "  " << i+size_ << "ssssss\n";}
   }

using index_t = int;
std::vector<double> mesh_val;
int size() {return size__;}
double operator()(int i) {return mesh_val[i];}
int index(int pos) {return pos;}
struct iterator {
       //iterator(){}
       iterator( Mesh<FERMI>& mesh_):mesh(mesh_), pos(0), size__(mesh_.size()){}
       void increment(){++pos;return;}
   Mesh<FERMI> & mesh;
  int pos,size__;
  bool is_end() { return (pos>size__-1); }
  int index() {return pos;}
  int r_index() {return pos-size__/2;}
  double val() { return this->mesh(pos);}
  void reset(){pos=0;}
  int size(){return size__;}
};
iterator begin(){
 return iterator(*this);
}

double T;
int size__;

};

template <>
struct Mesh<BOSON> {
using type = double;
Mesh(){}
Mesh(double T_,int size_): T(T_),mesh_val(size_){
    for (int i=0;i<mesh_val.size();++i) mesh_val[i] = pi*2.0*i*T;
   }
std::vector<double> mesh_val;
int size() {return mesh_val.size();}
double operator()(int i) {return mesh_val[i];}

struct iterator {
       iterator(Mesh<BOSON>& mesh_):mesh(mesh_), pos(-mesh_.size()), size__(mesh_.size()){}
       void increment(){++pos;return;}
  Mesh<BOSON> & mesh;
  int pos,size__;
  bool is_end() { return (pos>size__/2-1); }
  int index() {return size__/2+pos;}
  int val() { return mesh(pos);}
};

iterator begin(){
 return iterator(*this);
}


double T;
};

struct scalar{};
struct multi{};

template <typename MESH_TYPE, typename VALTYPE, typename DIM>
class GF{};

template <typename MESH_TYPE, typename VALTYPE>
class
GF<MESH_TYPE,VALTYPE, multi>{
public:
GF(MESH_TYPE mesh_, int N_):N(N_), mesh(mesh_), data(mesh_.size(), matrix<VALTYPE> (N_,N_, VALTYPE(0) ) ) {}
GF(){}
~GF(){}
using it_t = typename MESH_TYPE::iterator;
 matrix<VALTYPE>& operator()(typename MESH_TYPE::iterator& it){
  return data[ it.index()];
 }

 matrix<VALTYPE>& operator()(int i){
 return data[ this->mesh.size()/2 + i];
 }

int N;
MESH_TYPE mesh;
std::vector < matrix<VALTYPE> > data;
};

template <typename MESH_TYPE, typename VALTYPE>
class
GF<MESH_TYPE,VALTYPE, scalar>{
public:
GF(MESH_TYPE mesh_):mesh(mesh_), data(mesh_.size(), VALTYPE(0))  {}
GF(){}

using it_t = typename MESH_TYPE::iterator;
//GF(MESH_TYPE& mesh_, int& N_):N(N_), mesh(mesh_), data(mesh_.size(), matrix<VALTYPE>(N_,N_, VALTYPE(0)))  {}
~GF(){}
 VALTYPE& operator()(int i){
 return data[i];
 }

  VALTYPE& operator()(typename MESH_TYPE::iterator& it){
 return data[it.index()];
 }
int N;
MESH_TYPE mesh;
std::vector < VALTYPE > data;
};

template <typename A, typename MESH_TYPE, typename VALTYPE, typename DIM>
GF<MESH_TYPE,VALTYPE, DIM> operator*  (  A& a,  GF<MESH_TYPE,VALTYPE, DIM>& g1)
{
using GF_TYPE = GF<MESH_TYPE,VALTYPE, DIM>;
GF_TYPE tmp(g1);
for (typename GF_TYPE::it_t it(tmp.mesh); !it.is_end(); it.increment()){
   tmp(it) = a*g1(it) ;
}
return tmp;
}

  
template <typename A, typename MESH_TYPE, typename VALTYPE, typename DIM>
GF<MESH_TYPE,VALTYPE, DIM> operator*  ( A&& a,  GF<MESH_TYPE,VALTYPE, DIM>& g1)
{
using GF_TYPE = GF<MESH_TYPE,VALTYPE, DIM>;
GF_TYPE tmp(g1);
for (typename GF_TYPE::it_t it(tmp.mesh); !it.is_end(); it.increment()){
   tmp(it) = a*g1(it) ;
}
return tmp;
}
  

template <typename MESH_TYPE, typename VALTYPE, typename DIM>
GF<MESH_TYPE,VALTYPE, DIM> operator+  (  GF<MESH_TYPE,VALTYPE, DIM>&& g1,    GF<MESH_TYPE,VALTYPE, DIM>&& g2)
{
using GF_TYPE = GF<MESH_TYPE,VALTYPE, DIM>;
GF_TYPE tmp(g1);
for (typename GF_TYPE::it_t it(g1.mesh); !it.is_end(); it.increment()){
   tmp(it) = g1(it) + g2(it); 
   }
return tmp;
}


template <typename MESH_TYPE, typename VALTYPE>
class Kernel: public GF<MESH_TYPE,VALTYPE,scalar> {
public:
Kernel(){}
Kernel(MESH_TYPE mesh_, int nqstep_, std::string a2ffil_, double lowcut):GF<MESH_TYPE, VALTYPE, scalar>(mesh_), 
                                                          nqstep(nqstep_), a2f( nqstep, std::vector< double>(10,0.0))
                                                          ,w(nqstep,0.0)  
{
  std::ifstream a2ffil_stream(a2ffil_);
  std::vector<double> a2fw(nqstep,0.0);
  std::vector<double> w2(nqstep,0.0);
  std::vector<double> iw2(this->mesh.size(),0.0);
  std::ofstream of("a2f");
  for (int i=0; i< nqstep; ++i){
      a2ffil_stream >> w[i] >> a2f[i][0] >> a2f[i][1] >> a2f[i][2] ; //>> a2f[i][3] >> a2f[i][4] >> a2f[i][5] >> a2f[i][6]>>a2f[i][7]>> a2f[i][8]>>a2f[i][9];
      //of <<  w[i]<< " " << a2f[i][0] << "\n";
      w[i]*=0.001*1.0; // mev to ev
  }
  double weight = w[2]-w[1];
  for (int iw=0; iw < this->mesh.size(); ++iw)
      iw2[iw] = std::pow(this->mesh(iw),2.0);
  for ( int i=0; i< nqstep; ++i){
      a2fw[i] =  -2.0*w[i]*a2f[i][1]*weight;
/* minus sign is because of the form of non-int phonon propagator, i. e. there is 
         a term as
    2\omega/[(\imath v_n)^2- \omega^2]
  extra sign arises from \imath^2
*/
      w2[i] = std::pow(w[i],2.0);
  }
 std::ofstream oker("ker.dat"); 
  for (int iw=0; iw < this->mesh.size(); ++iw){
      for ( int i=0; i< nqstep; ++i) 
        if(w[i] > lowcut) this->data[iw] += a2fw[i]*std::pow(w2[i]+iw2[iw],-1.0);
      oker<<this->mesh(iw)<< " " << this->data[iw] << "\n";      
           
  } 
  std::cout<< lowcut << " Lambda(\Gamma(0)) =" <<  this->data[0] << "\n";
return;
}
private:
int nqstep;
std::vector< double >  w;
//std::vector< matrix<VALTYPE> > data;
std::vector< std::vector < double> > a2f;
};

struct first_diagram{};
struct second_diagram{};
  
template <typename DIAGORDER, typename DISPTYPE>
struct eliashberg;
  
template <typename T>
struct dispersion;

struct fromfile{};
struct square{};

//

struct Lmesh { 
Lmesh(){}
Lmesh(int size_,double min_, double max_):size__(size_),step((max_-min_)/size_), min__(min_),max__(max_){}
int size() {return size__;}
double operator() (int pos) {return min__ + step*pos;}
using index_t = int;
struct iterator {
       //iterator(){}
       iterator(Lmesh& mesh_):mesh(mesh_), pos(0), size__(mesh_.size()){}
       void increment(){++pos;return;}
   Lmesh & mesh;
    bool is_end() { return (pos>size__-1); }
  int index() {return pos;}
  int r_index() {return pos;} //-size__/2;}
  double val() { return mesh(pos);}
  void reset(){pos=0;}
  int size__,pos;
  int size() {return size__;}
};
iterator begin(){
 return iterator(*this);
}
double step,min__,max__;
int size__;
};       

struct Dmesh {
Dmesh(){}
Dmesh(int min_, int max_):size__(max_-min_),min__(min_), max__(max_){ }
int size() {return size__;}
using index_t = int;
double operator() (int pos) {return min__+pos;}

int index(int pos) {return pos;}
struct iterator {
       //iterator(){}
       iterator(Dmesh& mesh_):mesh(mesh_), pos(0), size__(mesh_.size()), min__(mesh_.min__){ }
       void increment(){++pos;return;}
  Dmesh & mesh;
    bool is_end() { return (pos>size__-1); }
   int index() {return (int) pos;}
  int r_index() {return pos+min__;}
  double val() { return mesh(pos);}
  void reset(){pos=0;}
  int size() {return size__;}
  int size__,pos,min__;
};
iterator begin(){
 return iterator(*this);
}

int min__,max__;
int size__;
};

template < int R = 0, typename ...T, typename ...F>
inline typename std::enable_if< R == sizeof...(T)-1, bool>::type increment__ (std::tuple<T...>& t, std::tuple<F...>& f)
{
using type = typename std::tuple_element<R, std::tuple<T...>>::type;
  std::get<R>(f).increment();
  
  if ( std::get<R>(f).is_end() ){
      return true;
  } else { std::get<R>(t) = std::get<R>(f).r_index();
           return false; }

}

template < int R = 0, typename ...T, typename ...F >
inline typename std::enable_if< R < sizeof...(T)-1  , bool>::type increment__ ( std::tuple<T...>& t, std::tuple<F...>& f)
{
  using type = typename std::tuple_element<R, std::tuple<T...>>::type;
  std::get<R>(f).increment();
  if ( std::get<R>(f).is_end() ){
      std::get<R>(f).reset();
      std::get<R>(t) = std::get<R>(f).r_index();
      increment__<R+1> ( t, f);
  } else { std::get<R>(t) = std::get<R>(f).r_index();return false;}

}

template < int R = 0, typename ...T, typename ...F>
inline typename std::enable_if< R == sizeof...(T)-1, void>::type do_ (std::tuple<T...>& t, std::tuple<F...>& f)
{
using type = typename std::tuple_element<R, std::tuple<T...>>::type;

   std::get<R>(t) = std::get<R>(f).r_index();

}

template < int R = 0, typename ...T, typename ...F >
inline typename std::enable_if< R < sizeof...(T)-1  , void>::type do_ ( std::tuple<T...>& t, std::tuple<F...>& f)
{
      std::get<R>(t) = std::get<R>(f).r_index();
      do_<R+1>(t,f);
}

struct index__{
template <typename T, typename U>
int  operator ()(T& t,U u){ return t.index()+t.size()*u;}
template <typename T>
int operator ()(T& t) {return t.index();}
};

struct size__{
template <typename T, typename U>
int  operator ()(T& t,U u){ return t.size()*u;}
template <typename T>
int operator ()(T& t) {return t.size();}
};


template < int R = 0, typename T, typename ...F>
inline typename std::enable_if< R == sizeof...(F)-1, int>::type do__ ( T& t, std::tuple<F...>& f)
{ 
return t(std::get<R>(f)); 
}

template < int R = 0, typename T, typename ...F >
inline typename std::enable_if< R < sizeof...(F)-1  , int>::type do__ ( T& t, std::tuple<F...>& f)
{
  return t(std::get<R>(f), do__<R+1>(t,f));
}


template < int R, typename T>
int p(T& t) {return std::get<R>(t.get_r_index_tuple());}

template < int R, typename T>
double pv(T& t) {return std::get<R> (t.get_it_tuple()).val();}

template < int R, typename T>
int idx__(T& t) {return std::get<R> (t.get_it_tuple()).index();}

template < int R, typename T>
typename std::tuple_element<R,typename T::it_tuple_t>::type& it__(T& t) {return std::get<R> (t.get_it_tuple());}

template <typename... Ms> struct mesh_product{
  mesh_product(){}
   mesh_product(Ms&... t_):m_tuple(std::tuple<Ms...>(t_...)), it_tuple( it_tuple_t(t_.begin()...)) {}
   mesh_product(Ms&&... t_):m_tuple(std::tuple<Ms&...>(t_...)), it_tuple( it_tuple_t(t_.begin()...)) {}
  using m_tuple_t = std::tuple<Ms...>;
  using it_tuple_t = std::tuple<typename Ms::iterator...>;
  m_tuple_t m_tuple;
  it_tuple_t it_tuple;
  static constexpr size_t tuple_size = sizeof...(Ms);
  m_tuple_t get_m_tuple() {return m_tuple;} 
  size__ size_temp; 
  int size() {return do__(size_temp,it_tuple);}  
  it_tuple_t begin() { return it_tuple;} 
// ***************************************
struct iterator {
       iterator ( mesh_product<Ms...>& mp_):mp(mp_),it_tuple(mp_.begin())
       {
         do_(r_index_tuple, it_tuple);
       }
       mesh_product<Ms...>& mp;
       index__ index_tmp;
       size__ size_temp;
       using it_tuple_t = std::tuple<typename Ms::iterator...>;
       std::tuple<typename Ms::index_t...> r_index_tuple; 
       std::tuple<typename Ms::iterator...> it_tuple;
       bool increment () { return increment__(r_index_tuple,it_tuple);}
       std::tuple<typename Ms::index_t...>& get_r_index_tuple () {return r_index_tuple;}
       std::tuple<typename Ms::iterator...>& get_it_tuple() { return it_tuple;}
       bool is_end() { return std::get<tuple_size-1>(it_tuple).is_end();}
       int index() {return do__(index_tmp,it_tuple); }
       int size() {return do__(size_temp,it_tuple);}
       template<int R>
       int   p(){return std::get<R>(r_index_tuple);}
};

};

template <>
struct dispersion<square> {
dispersion(){}
dispersion(double t_, int size_):size__(size_), t(t_), step(2.0*pi/size_) {}
 using it_t = typename mesh_product < Lmesh, Lmesh > :: iterator ;
int size(){ return size__; }
 double operator ()(int kx, int ky)  { return  -2.0*t*( cos(step*kx) + cos(step*ky) + 0.0*t*sin(step*kx)*sin(step*ky)  ); }

 double operator () (it_t& it ){ return  -2.0*t*( cos(step*p<0>(it)) + cos(step*p<1>(it)) );}

double dk() { return step;}
private:
int size__;
double t;
double step;
int nbands;
};

template <>
struct dispersion<fromfile> {
dispersion(){}

using it_t = typename mesh_product < Lmesh, Lmesh > :: iterator ;
dispersion( std::string bandfile_, int nbands_, double Ef_, double sigma_, int Nk_):e(0, std::vector<double>(nbands_,0.0)), klist (0, std::vector<double> (2,0.0)), dos(nbands_, 0), Ef(Ef_),nbands(nbands_), sigma(sigma_), Nk__(Nk_){
   std::ifstream fi(bandfile_ );
   std::string s;
   int npoints;
   
   std::vector<double> ktmp (2,0.0);
   std::vector<double> etmp (8,0.0);
   double tmp = 10, Evbm(-1000.0);
   npoints=Nk_*Nk_;
//   fi >> s >> s >> nbands >> s >> s >> npoints ;
// std::cout << "Nbands = " << nbands << "npoits = " << npoints << "\n";
   //size__ = sqrt(npoints);
   for ( int i=0 ; i<npoints; ++i){
   fi >> ktmp[0] >> ktmp[1] >> tmp;
// std::cout << ktmp[0] << " "<< ktmp[1] << "\n";
   for ( int j=0 ; j<8; ++j) {
       fi >> etmp[j];
       if ( j < 5 && etmp[j] > Evbm ) Evbm = etmp[j];
   }

   klist.push_back(ktmp);
   e.push_back(etmp);
   }
 {
 double tmp_dos(0);
 using kmesh2_t =  mesh_product < Lmesh, Lmesh > ;
 Lmesh kmesh(Nk_,0,2*pi);
 Delta delta;
 kmesh2_t kmesh2(kmesh,kmesh);
 double weight = std::pow(kmesh.size(),-2.0);
 std::cout << " Evbm= " << Evbm << " DeltaEf = " << Ef-Evbm << "\n";
 for (int i=0; i< nbands; ++i){
    for ( it_t  kit(kmesh2) ; !kit.is_end(); kit.increment() ){
        tmp_dos+= weight*delta ((e[kit.index()][i]-Ef)/sigma )/sigma;
    }
    std::cout << " Dos of band " << i+1 << " = " << tmp_dos << "\n";
    dos[i] = tmp_dos;
    tmp_dos=0;
 }
// calculating DOS for a range of energies

std::vector<double> DOS(5000,0.0);
double wstep= (2.0-(-18.0))/5000.0, w(0.0); 
for ( int iw=0; iw<5000; iw++){
    w = -18.0+iw*wstep;
      for (int i=0; i< nbands; ++i){
          for ( it_t  kit(kmesh2) ; !kit.is_end(); kit.increment() ){
          DOS[iw]+= weight*delta ((e[kit.index()][i]-w)/sigma )/sigma;
          }
      }
    }
std::string dosfil = std::string("DOS_")+std::string("Nk")+std::to_string(Nk_) + std::string("sigma") +std::to_string( sigma ) + ".dat";
ofstream of(dosfil);
for ( int iw=0; iw<5000; iw++){
    w = -18.0+iw*wstep;
    of<< w << " " << DOS[iw] << "\n";
}

}
//
}
int size(){ return Nk__; }
std::vector<double> operator () ( it_t& it ){ return  e[it.index()]; }

//private:
int Nk__;
double Ef,sigma;
double step;
int nbands,band;
std::vector< std::vector <double> > klist,e;
std::vector <double> dos;
};



 
template <typename VALTYPE, typename DIM, typename ...MESH_TYPES>
class GF_NONLOC{};
 
template <typename VALTYPE>
class
GF_NONLOC< Mesh<FERMI>, Dmesh, Lmesh ,VALTYPE, multi>{
public:
GF_NONLOC( Mesh<FERMI>&  fmesh_, Dmesh& dmesh_, Lmesh& lmesh_, int N_):dim(N_), 
        rf_mesh( mesh_product<Dmesh, Dmesh, Mesh<FERMI> >(dmesh_,dmesh_,fmesh_) ),  
         kf_mesh( mesh_product<Lmesh, Lmesh, Mesh<FERMI> >(lmesh_,lmesh_,fmesh_) ), 
          rmesh(dmesh_), kmesh(lmesh_), fmesh(fmesh_), 
           r2_mesh(dmesh_,dmesh_), k2_mesh(lmesh_,lmesh_),
            data(fmesh_.size()*lmesh_.size()*lmesh_.size(), matrix<VALTYPE> (N_,N_, VALTYPE(0))), 
             data_r( fmesh_.size()*dmesh_.size()*dmesh_.size(), matrix<VALTYPE> (N_,N_, VALTYPE(0))), 
              Nw2(fmesh_.size()/2), Nk2(dmesh_.size()/2),Nr2(lmesh_.size()/2)
                  

{
zero_r_point =  dmesh_.size()*dmesh_.size()/2 + dmesh_.size();
}
~GF_NONLOC()
{

}
 
using  rf_mesh_t = mesh_product < Dmesh, Dmesh, Mesh<FERMI> >;
using  kf_mesh_t = mesh_product < Lmesh, Lmesh, Mesh<FERMI> >;
using  r2_mesh_t = mesh_product < Dmesh, Dmesh > ;
using  k2_mesh_t = mesh_product < Lmesh, Lmesh > ;
  matrix<VALTYPE>& Gr(int iw, int rx, int ry){
                return data_r[Nk2+ry+rmesh.size()*( rx+Nk2+rmesh.size()*(Nw2+iw))];
    }

  matrix<VALTYPE>& Gk(int iw, int kx, int ky){ // we suppose that k-mesh is defined in positive region of first brilium zone
                return data[ky+kmesh.size()*( kx+kmesh.size()*(Nw2+iw)) ];
    }


kf_mesh_t& get_kf_mesh() {return kf_mesh;}
 
matrix<VALTYPE>& Gr( typename rf_mesh_t::iterator& rit){
  return data_r[rit.index()];
}
 
matrix<VALTYPE>& Gk( typename kf_mesh_t::iterator& kit){
  return data[kit.index()];
}

matrix<VALTYPE>& Gk( typename Mesh<FERMI>::iterator fit, typename k2_mesh_t::iterator& kit){
  return data[fit.index()*kit.size() + kit.index()];
}

matrix<VALTYPE>& operator () ( typename kf_mesh_t::iterator& kit){
  return data[kit.index()];
}

matrix<VALTYPE>& operator () ( typename rf_mesh_t::iterator& rit){
  return data_r[rit.index()];
}


void fourier()
{
    Dmesh dm(0,dim);
    mesh_product < Dmesh, Dmesh > d(dm,dm); 

    std::complex<double>* in = new std::complex<double>[k2_mesh.size()];
    std::complex<double>* out = new std::complex<double>[r2_mesh.size()];
    fftw_plan fftp;
    std::ofstream of("Gr_sp.dat");
    double t = clock();
    std::complex<double> fac( std::exp(-imath*pi*(double) kmesh.size()/2.0 )   );
    for ( mesh_product < Dmesh, Dmesh >::iterator dit(d); !dit.is_end(); dit.increment() ) {
        for (Mesh<FERMI>::iterator fit(fmesh); !fit.is_end(); fit.increment()) {
            for (mesh_product < Lmesh, Lmesh >::iterator kit(k2_mesh); !kit.is_end(); kit.increment() )  
                in[kit.index()] = data[ fit.index() * kit.size() + kit.index()][p<0>(dit)][p<1>(dit)]*pow(kit.size(), -1.0)*std::exp(-imath*pi*(double)(idx__<0>(kit)+idx__<1>(kit)))*fac;
        fftp = fftw_plan_dft_2d(kmesh.size(), kmesh.size(), reinterpret_cast<fftw_complex*>( in) , reinterpret_cast<fftw_complex*> (out), FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(fftp); 
        for (mesh_product < Dmesh, Dmesh >::iterator rit(r2_mesh); !rit.is_end(); rit.increment() )  
            data_r[ fit.index() * rit.size() + rit.index()][p<0>(dit)][p<1>(dit)] = out[rit.index()] ; //* std::exp(-imath*pi*(double)(idx__<0>(rit)+idx__<1>(rit)));                        
    }
    fftw_destroy_plan(fftp);
  }
  t  = clock() - t;
  printf ("Fourier transform is done in (%f seconds).\n",((float)t)/CLOCKS_PER_SEC);

std::cout << r2_mesh.size() << "\n";
  for (Mesh<FERMI>::iterator fit(fmesh); !fit.is_end(); fit.increment()) {
      for (mesh_product < Dmesh, Dmesh >::iterator rit(r2_mesh); !rit.is_end(); rit.increment() )      
          of << fit.val() << " " << p<0>(rit) << " " << p<1>(rit) << " " << data_r[ fit.index()*rit.size() + rit.index() ][0][0].real() << " " << data_r[ fit.index()*rit.size() + rit.index() ][0][0].imag()<< " " << data_r[ fit.index()*rit.size() + rit.index() ][0][1].real() <<"  "<< data_r[ fit.index()*rit.size() + rit.index() ][0][1].imag() << "\n";
              
  }

ofstream of2("G_loc2.dat");
 int nk = rmesh.size();
 for (Mesh<FERMI>::iterator fit(fmesh); !fit.is_end(); fit.increment()) {
     of2 << fit.val() << " " << this->Gr(fit.r_index(), 0 , 0)[0][0] << " " << this->Gr(fit.r_index(), 0 , 0)[1][0] << "\n" ;
 } 


 delete[] in;
 delete[] out;
}
 
matrix<VALTYPE> G_loc(Mesh<FERMI>::iterator &fit){
return data_r[ fit.index() * r2_mesh.size() + zero_r_point];
} 

void fourier_back()
{
    mesh_product < Dmesh, Dmesh > rmesh2(rmesh,rmesh);
    mesh_product < Lmesh, Lmesh > kmesh2(kmesh,kmesh);
    Dmesh dm(0,dim);
    mesh_product < Dmesh, Dmesh > d(dm,dm);
    std::complex<double> fac( std::exp(imath*pi*(double) kmesh.size()/2.0 )   );
    std::complex<double>* in = new std::complex<double>[kmesh2.size()];
    std::complex<double>* out = new std::complex<double>[rmesh2.size()];
    fftw_plan fftp;
    double t = clock();
  for ( mesh_product < Dmesh, Dmesh >::iterator dit(d); !dit.is_end(); dit.increment() ) {
    for (Mesh<FERMI>::iterator fit(fmesh); !fit.is_end(); fit.increment()) {
      for (mesh_product < Dmesh, Dmesh >::iterator rit(rmesh2); !rit.is_end(); rit.increment() )
          in[rit.index()] = data_r[ fit.index() * rit.size() + rit.index()][p<0>(dit)][p<1>(dit)] ;

          fftp = fftw_plan_dft_2d(kmesh.size(), kmesh.size(), reinterpret_cast<fftw_complex*>( in), reinterpret_cast<fftw_complex*> (out), FFTW_BACKWARD, FFTW_ESTIMATE);
      fftw_execute(fftp);

      for (mesh_product < Lmesh, Lmesh >::iterator kit(kmesh2); !kit.is_end(); kit.increment() )
          data[ fit.index() * kit.size() + kit.index()][p<0>(dit)][p<1>(dit)] = out[kit.index()]*std::exp(imath*pi*(double)(idx__<0>(kit)+idx__<1>(kit)))*fac;
    }
    fftw_destroy_plan(fftp);
  }
  t  = clock() - t;
  printf ("Fourier back in done in (%f seconds).\n",((float)t)/CLOCKS_PER_SEC);
 
 delete[] in;
 delete[] out;
}

int dim;
rf_mesh_t rf_mesh;
kf_mesh_t kf_mesh;
Dmesh rmesh;
Lmesh kmesh;
k2_mesh_t k2_mesh;
r2_mesh_t r2_mesh;
Mesh < FERMI > fmesh;
std::vector < matrix<VALTYPE> > data,data_r;
int Nr2, Nk2, Nw2, zero_r_point;
};
   
template <typename GNONLOC0, typename SLOC , typename GNONLOC>
void Dyson (GNONLOC0& G0, SLOC& Sigma, GNONLOC& G){
  for ( typename GNONLOC::kf_mesh_t::iterator it(G0.get_kf_mesh()); !it.is_end(); it.increment() ){
      G(it) = Inv( G0(it) ) - Sigma(it__<2>(it));
      G(it) = Inv ( G(it) );
     // std::cout << G(it)[0][1] << " " << G(it)[1][0] <<" "<< Sigma(it__<2>(it))[0][1] << " "<< Sigma(it__<2>(it))[1][0]<<  "\n ";
  }
return;
} 
template <typename VALTYPE, typename DISPTYPE>
class G0_:public GF_NONLOC< Mesh<FERMI>, Dmesh, Lmesh ,VALTYPE, multi> {
public:
G0_(Mesh<FERMI>&  fmesh_, Dmesh& dmesh_, Lmesh& lmesh_, int N_, DISPTYPE dispersion_, double Ef_, int band_):GF_NONLOC< Mesh<FERMI>, Dmesh, Lmesh ,VALTYPE, multi>(fmesh_,dmesh_,lmesh_,N_), ek ( dispersion_)
{ 
std::ofstream of("Gk_sp.dat");
 for (typename Mesh<FERMI>::iterator fit(this->fmesh); !fit.is_end(); fit.increment() )
     for ( typename GF_NONLOC< Mesh<FERMI>, Dmesh, Lmesh ,VALTYPE, multi>::k2_mesh_t::iterator kit(this->k2_mesh) ; !kit.is_end(); kit.increment() ){ 
//      cout << std::get<2>(it.get_it_tuple()).val() << " "<< std::get<2>(it.get_it_tuple()).index()  << "\n";
      this->data[fit.index()*kit.size()+kit.index()][0][0] =  pow( imath*fit.val()+ Ef_ - ek( kit)[band_], -1.0 );
      this->data[fit.index()*kit.size()+kit.index()][1][1] =  pow( imath*fit.val()- Ef_ + ek( kit)[band_], -1.0 );
    }
}

 DISPTYPE ek;
};  


void Lambda(std::string a2ffil, int nqstep, double sctemp, int Niw, double muc,double phcut_, double& w2, double& wn, double& mucn)
{
  std::vector< std::vector < double > > a2f( nqstep, std::vector< double >(5,0.0));
  std::vector< double >  w(nqstep,0.0);
  std::ifstream a2ffil_stream(a2ffil);
  w2 = 0.0; wn = 0.0; mucn = 0.0;

  for (int i=0; i< nqstep; ++i){
      a2ffil_stream >> w[i] >> a2f[i][0] >> a2f[i][1] >> a2f[i][2];// >>  a2f[i][3] >>  a2f[i][4];
      w[i]*=0.001*1.0; // mev to ev
  }
  double weight( w[2]-w[1] ),lambda(0);

  for (int i=0; i< nqstep; ++i){
    if (w[i] > phcut_ ) {
      lambda += 2.0*a2f[i][0]/w[i];
      w2 += 2.0*a2f[i][0] * w[i] ;
      }
  }

  lambda*=weight;
  w2 *= (weight / lambda);
  wn  = sctemp*pi*( Niw + 1.0 );

  w2 = std::pow(w2,0.5);

  mucn = muc/( 1.0 + muc*std::log(w2/wn) );
  std::cout<<"mucn ="<<mucn << " lambda=  " << lambda << "w2 ="<<w2<<" wn = "<<wn<< "\n ";

 return;
}



template <typename DISPSRC>
struct eliashberg<second_diagram,DISPSRC> {         
public:
eliashberg(){}
eliashberg(double T_, double muc_, int Niw_, int Nk_, double initial_gap_, bool second_diag_, string a2ffil_, dispersion<DISPSRC> dispersion_, double phcut_, std::string basename_, int Nr_): Niw(Niw_), T(T_), phcut(phcut_),muc(muc_),  phi(Mesh<FERMI>(T_, 3*Niw_ )), chi(Mesh<FERMI>(T_, 3*Niw_ )), 
                                                              zi(Mesh<FERMI>(T_, 3*Niw_ )), S(Mesh<FERMI>(T_, 3*Niw_ ), 2 ),
                                                              Ker( Mesh<BOSON>(T_, 4.0*(Niw_+2)),2000, a2ffil_, phcut_ ),
                                                              s0(2,2,czero),s1(2,2,czero),s2(2,2,czero),s3(2,2,czero),
                                                              maxiter(500), fmesh(T_, Niw_), kmesh(Nk_,0,2*pi),
                                                              rmesh(-Nk_/2,Nk_/2), ek(dispersion_), second_diag(second_diag_),
                                                              basename(basename_), Nr(Nr_)
{
s0[0][0]=cone;s0[1][1]=cone;
s1[0][1]=cone;s1[1][0]=cone;
s2[0][1]=imath;s2[1][0]=-imath;
s3[0][0]=cone;s3[1][1]=-cone;

bool flag(true);


for ( Mesh<FERMI>::iterator w(fmesh); !w.is_end(); w.increment() ) {
    S(w.r_index())[0][1] = initial_gap_; S(w.r_index())[1][0] = initial_gap_;
    S(w.r_index())[0][0] = czero; S(w.r_index())[1][1] = czero;
}

Lambda(a2ffil_, 2000, T_, Niw_ , muc_, phcut_, w2, wn, mucn);
}

~eliashberg(){} 
void eval_Sigma(){
     std::complex<double> tmp(0,0);
     std::complex<double> zeta(czero),zi_tmp(czero), phi_tmp(czero);
     Mesh<FERMI> fmesh(T, Niw);
     GF<Mesh<FERMI>, std::complex<double>, multi > S1(fmesh,2),S2(fmesh,2);
     Mesh<FERMI> fmesh3(T, 3*Niw );     
     GF_NONLOC< Mesh<FERMI>, Dmesh, Lmesh , std::complex<double>, multi> g_(fmesh3,rmesh,kmesh,2), s_(fmesh,rmesh,kmesh,2);
     G0_<std::complex<double>, dispersion<DISPSRC> > g0(fmesh3,rmesh,kmesh,2, ek , ek.Ef,4);
  
     matrix < std::complex<double> > A(2,2,fczero),B(2,2,fczero),C(2,2,fczero);
     std::complex<double> temp3(czero), temp4(czero); //, kern(fczero);
     double kern;
// self energy is local
     Dyson(g0, S, g_);
// fourier transform to real space
     g_.fourier();
     double t = clock();
// debug for the first order k-summation comparing to FFT
#ifdef KSUM
GF<Mesh<FERMI>, std::complex<double>, multi > G_local(Mesh<FERMI>(T, 3*Niw ), 2 );
{
mesh_product< Lmesh, Lmesh> k2mesh (kmesh, kmesh);
for ( Mesh<FERMI>::iterator w(fmesh3); !w.is_end(); w.increment() ){
    for ( mesh_product< Lmesh, Lmesh>::iterator kit(k2mesh); !kit.is_end(); kit.increment())
        G_local( w.r_index() )+= pow(k2mesh.size(),-1.0)*g_.Gk(w.r_index(), p<0>(kit), p<1>(kit));
    }
}

for ( Mesh<FERMI>::iterator w(fmesh3); !w.is_end(); w.increment() ){
    std::cout<<w.val()<< "  direct G_00 = " << G_local( w.r_index() )[0][0] << "fft G_00 = " << g_.Gr( w.r_index(), 0 ,0)[0][0] << "\n";
    std::cout<<w.val() << "  direct G_01 = " << G_local( w.r_index() )[0][1] << "fft G_01 = " << g_.Gr( w.r_index(), 0 ,0)[0][1] << "\n";
    std::cout<<w.val() << "  direct G_10 = " << G_local( w.r_index() )[1][0] << "fft G_10 = " << g_.Gr( w.r_index(), 0 ,0)[1][0] << "\n";
    std::cout<<w.val() << "  direct G_11 = " << G_local( w.r_index() )[1][1] << "fft G_11 = " << g_.Gr( w.r_index(), 0 ,0)[1][1] << "\n";
}
#endif
  
// first diagram evaluation
     mesh_product < Mesh<FERMI>, Mesh<FERMI> > m_w2(fmesh,fmesh);
     for (mesh_product < Mesh<FERMI>, Mesh<FERMI> >::iterator w(m_w2); !w.is_end(); w.increment() ){
            S1( p<0>(w) )  +=  pow(ek.dos[4],-1.0) *  Ker( abs( p<0>(w) - p<1>(w)) ) * g_.Gr( p<1>(w), 0 ,0); //G_local is G_{r=0}(iw) 
            S1( p<0>(w) )[0][1] +=  pow(ek.dos[4],-1.0) * mucn * g_.Gr( p<1>(w), 0, 0)[0][1] ; // G_local is G_{r=0}(iw)
            S1( p<0>(w) )[1][0] +=  pow(ek.dos[4],-1.0) * mucn * g_.Gr( p<1>(w), 0, 0)[1][0] ;
     }
     for ( Mesh<FERMI>::iterator w(fmesh); !w.is_end(); w.increment() ){
         S1 ( w.r_index() ) = (-T) * s3 * ( S1( w.r_index() ) * s3 );
         S ( w.r_index() )=  S1( w.r_index() );
         
     if ( w.r_index() == 0)
        std::cout << "Delta(w0) extracted from the first diagramm : "<< w.val() << " " << (( S(w.r_index())[0][1] + S(w.r_index())[1][0] ).real() * 0.5)/(1.0 - S(w.r_index())[0][0].imag()/w.val()) << "\n";
     }
     std::cout << " lllllllllllllll " <<  g_.Gr(0, 0 , 0)[0][0] << "\n";
     
     t = clock() - t;
     printf ("Frequency summation for first diagram (%f seconds).\n",((float)t)/CLOCKS_PER_SEC); 
     
     if ( second_diag == true ) {
        t = clock(); 
        mesh_product<Mesh<FERMI> , Mesh<FERMI>, Mesh<FERMI> > iw_mesh(fmesh, fmesh, fmesh);
        mesh_product< Dmesh, Dmesh> r2mesh (rmesh, rmesh);
        mesh_product< Lmesh, Lmesh> k2mesh (kmesh, kmesh);
// second diagramm evaluation
     double dos4pm2 = std::pow(ek.dos[4], -2.0);
     if ( 1 == 0 ){
        for (mesh_product<Mesh<FERMI> , Mesh<FERMI>, Mesh<FERMI> >::iterator w(iw_mesh); !w.is_end(); w.increment() ) {
           for ( mesh_product< Dmesh, Dmesh>::iterator x(r2mesh); !x.is_end(); x.increment() )
           s_.Gr(p<0>(w), p<0>(x) , p<1>(x)) += dos4pm2 * Ker( abs(p<0>(w)-p<1>(w) ) ) * Ker( abs( p<1>(w) - p<2>(w) ) ) * 
                                                                            g_.Gr( p<1>(w), p<0>(x) , p<1>(x) ) * 
                                                                  ( s3 * g_.Gr(p<2>(w), -p<0>(x) , -p<1>(x) ) ) * 
                                               (  s3 * g_.Gr(p<0>(w) + p<2>(w) - p<1>(w), p<0>(x) , p<1>(x) ) ) ;
  
        }
  
      } else {
        Dmesh rmesh_cut( -Nr, Nr );
        mesh_product< Dmesh, Dmesh> rmesh2_cut(rmesh_cut,rmesh_cut);
        
        for (mesh_product<Mesh<FERMI> , Mesh<FERMI>, Mesh<FERMI> >::iterator w(iw_mesh); !w.is_end(); w.increment() ) {
           for ( mesh_product< Dmesh, Dmesh>::iterator x(rmesh2_cut); !x.is_end(); x.increment() ){
              for (int i=0; i<2; ++i)
                 for ( int j =0 ; j<2; ++j){
                     A[i][j] = g_.Gr( p<1>(w), p<0>(x) , p<1>(x) )[i][j]; 
                     B[i][j] = g_.Gr(p<2>(w), -p<0>(x) , -p<1>(x) )[i][j] ; 
                     C[i][j] = g_.Gr(p<0>(w) + p<2>(w) - p<1>(w), p<0>(x) , p<1>(x) )[i][j] ;
                 }
                temp3 =  A[0][0] * B[0][0] - B[1][0] * A[0][1];
                temp4 = -B[0][1] * A[0][0] + B[1][1] * A[0][1];
                kern =  ( Ker( abs(p<0>(w)-p<1>(w) ) ) * Ker( abs( p<1>(w) - p<2>(w) ) ) ) * dos4pm2;
                s_.Gr(p<0>(w), p<0>(x) , p<1>(x))[0][0] += (C[0][0] * temp3 + C[1][0] * temp4 ) * kern;
                s_.Gr(p<0>(w), p<0>(x) , p<1>(x))[0][1] += (C[0][1] * temp3 + C[1][1] * temp4 ) * kern;
              
          }
        }
      }
        t  = clock() - t;
        printf ("Frequency summation for second diagram (%f seconds).\n",((float)t)/CLOCKS_PER_SEC);

        s_.fourier_back();
        
        Delta delta;
        double sigma_ = ek.sigma;
        double Tp2 = pow(T,2.0); 
        for (Mesh<FERMI>::iterator f(fmesh); !f.is_end(); f.increment() ){
           for ( mesh_product< Lmesh, Lmesh>::iterator k(k2mesh); !k.is_end(); k.increment() ) {
               S2( f.r_index() ) += delta( (ek( k )[4]-ek.Ef)/sigma_ )/sigma_/k2mesh.size()/ek.dos[4] * s_.Gk( f.r_index(), p<0>(k), p<1>(k) );
     //            if ( p<0>(k) == 0 && p<1>(k) == 0 )
     //                std::cout<< f.val() << " " << p<0>(k) << " " << p<1>(k)<<" " << s_.Gk(f.r_index() , p<0>(k) , p<1>(k))[0][0] << " " << s_.Gk(f.r_index(), p<0>(k) , p<1>(k))[0][1] << "\n";
           }
           S2( f.r_index() )[1][1] = -std::conj( S2( f.r_index() )[0][0] ); S2( f.r_index() )[1][0] = S2( f.r_index() )[0][1];
            S2( f.r_index() ) = s3 * ( Tp2 * S2( f.r_index() ) * s3);
           S( f.r_index() ) += S2( f.r_index() ); //s3 * ( Tp2 * S2( f.r_index() ) * s3);
           }
    }

   ofstream ioker(std::string("ISigma_kavg_") + basename + std::to_string(iter));
   ofstream roker(std::string("RSigma_kavg_") + basename + std::to_string(iter));

   ofstream ioker1(std::string("ISigma1_kavg_") + basename + std::to_string(iter));
   ofstream roker1(std::string("RSigma1_kavg_") + basename + std::to_string(iter));

   ofstream ioker2(std::string("ISigma2_kavg_") + basename + std::to_string(iter));
   ofstream roker2(std::string("RSigma2_kavg_") + basename + std::to_string(iter));


   for (Mesh<FERMI>::iterator f(fmesh); !f.is_end(); f.increment() ){
       ioker << f.val() << " " << S( f.r_index() )[0][0].imag() << " " << S( f.r_index() )[0][1].imag() << " " << S( f.r_index() )[1][0].imag() <<  " " << S( f.r_index() )[1][1].imag() << "\n";
           roker << f.val() << " " << S( f.r_index() )[0][0].real() << " " << S( f.r_index() )[0][1].real() << " " << S( f.r_index() )[1][0].real() <<  " " << S( f.r_index() )[1][1].real() << "\n";
 
   }

   for (Mesh<FERMI>::iterator f(fmesh); !f.is_end(); f.increment() ){
       ioker1<< f.val() << " " << S1( f.r_index() )[0][0].imag() << " " << S1( f.r_index() )[0][1].imag() << " " << S1( f.r_index() )[1][0].imag() <<  " " << S1( f.r_index() )[1][1].imag() << "\n";
       roker1<< f.val() << " " << S1( f.r_index() )[0][0].real() << " " << S1( f.r_index() )[0][1].real() << " " << S1( f.r_index() )[1][0].real() <<  " " << S1( f.r_index() )[1][1].real() << "\n";

   }

   for (Mesh<FERMI>::iterator f(fmesh); !f.is_end(); f.increment() ){
       ioker2<< f.val() << " " << S2( f.r_index() )[0][0].imag() << " " << S2( f.r_index() )[0][1].imag() << " " << S2( f.r_index() )[1][0].imag() <<  " " << S2( f.r_index() )[1][1].imag() << "\n";
       roker2<< f.val() << " " << S2( f.r_index() )[0][0].real() << " " << S2( f.r_index() )[0][1].real() << " " << S2( f.r_index() )[1][0].real() <<  " " << S2( f.r_index() )[1][1].real() << "\n";
   }

   ioker.close();
   roker.close();

   ioker1.close();
   roker1.close();

   ioker2.close();
   roker2.close();


return;
}

       int back_to_bz( int  k){
          if ( k< 0)
            if ( k+(-k/kmesh.size())*kmesh.size() < 0 )
               return k+(-k/kmesh.size())*kmesh.size() + kmesh.size();
            else 
               return 0;
          else 
            return k;
       }
void eval_Sigma_direct(){
     ofstream oker("sigma_k");
     std::complex<double> tmp(0,0);
     std::complex<double> zeta(czero),zi_tmp(czero), phi_tmp(czero);
     Mesh<FERMI> fmesh(T, Niw);
     GF<Mesh<FERMI>, std::complex<double>, multi > S1(fmesh,2),S2(fmesh,2);
     Mesh<FERMI> fmesh3(T, 3*Niw );     
     GF_NONLOC< Mesh<FERMI>, Dmesh, Lmesh , std::complex<double>, multi> g_(fmesh3,rmesh,kmesh,2), s_(fmesh,rmesh,kmesh,2);
     G0_<std::complex<double>, dispersion<DISPSRC> > g0(fmesh3,rmesh,kmesh,2, ek , ek.Ef,4);
  
     matrix < std::complex<double> > A(2,2,fczero),B(2,2,fczero),C(2,2,fczero);
     std::complex<double> temp3(czero), temp4(czero); //, kern(fczero);
     double kern;
// self energy is local
     cout << S(0)[0][1] << "\n"; 
     Dyson(g0, S, g_);

// fourier transform to real space to evaluate local G
     g_.fourier();
     double t = clock();
// first diagram evaluation
     mesh_product < Mesh<FERMI>, Mesh<FERMI> > m_w2(fmesh,fmesh);
     for (mesh_product < Mesh<FERMI>, Mesh<FERMI> >::iterator w(m_w2); !w.is_end(); w.increment() ){
            S1( p<0>(w) )  +=  pow(ek.dos[4],-1.0) *  Ker( abs( p<0>(w) - p<1>(w)) ) * g_.Gr( p<1>(w), 0 ,0); //G_local is G_{r=0}(iw) 
            S1( p<0>(w) )[0][1] += pow(ek.dos[4],-1.0) * mucn * g_.Gr( p<1>(w), 0, 0)[0][1] ; // G_local is G_{r=0}(iw)
            S1( p<0>(w) )[1][0] += pow(ek.dos[4],-1.0) * mucn * g_.Gr( p<1>(w), 0, 0)[1][0] ;
     }
       
     for ( Mesh<FERMI>::iterator w(fmesh); !w.is_end(); w.increment() ){
         S ( w.r_index() ) = s3 * ( S1( w.r_index() ) * s3 );

         S ( w.r_index() )= (-T) * S( w.r_index() );

     if ( w.r_index() == 0)
        std::cout << "Delta(w0) extracted from the first diagramm : "<< w.val() << " " << (( S(w.r_index())[0][1] + S(w.r_index())[1][0] ).real() * 0.5)/(1.0 - S(w.r_index())[0][0].imag()/w.val()) << "\n";
     }
     oker.close();
     t = clock() - t;
     printf ("Frequency summation for first diagram (%f seconds).\n",((float)t)/CLOCKS_PER_SEC); 
     
     if ( second_diag == true ) {
        t = clock(); 
        mesh_product<Mesh<FERMI> , Mesh<FERMI>, Mesh<FERMI> > iw_mesh(fmesh, fmesh, fmesh);
        mesh_product< Dmesh, Dmesh> r2mesh (rmesh, rmesh);
        mesh_product< Lmesh, Lmesh> k2mesh (kmesh, kmesh), q2mesh (kmesh, kmesh), qprim2mesh (kmesh, kmesh);
 
     if ( 1 == 0 ){
         for ( mesh_product< Lmesh, Lmesh>::iterator k(k2mesh); !k.is_end(); k.increment() ){
           for (mesh_product<Mesh<FERMI> , Mesh<FERMI>, Mesh<FERMI> >::iterator w(iw_mesh); !w.is_end(); w.increment() ) {
              for ( mesh_product< Lmesh, Lmesh>::iterator q(q2mesh); !q.is_end(); q.increment() )
                  for ( mesh_product< Lmesh, Lmesh>::iterator qpr(k2mesh); !qpr.is_end(); qpr.increment() )
           s_.Gk(p<0>(w), p<0>(k) , p<1>(k)) += std::pow(k2mesh.size(),-2.0)* Ker( abs(p<0>(w)-p<1>(w) ) ) * Ker( abs( p<1>(w) - p<2>(w) ) ) * 
                                                                            g_.Gk( p<1>(w), back_to_bz(p<0>(k)-p<0>(q)) ,back_to_bz( p<1>(k)-p<1>(q)) ) * 
                                                                  ( s3 * g_.Gk(p<2>(w), back_to_bz(p<0>(k)-p<0>(q)-p<0>(qpr)), back_to_bz(p<1>(k)-p<1>(q)-p<1>(qpr) ) ) ) * 
                                               (  s3 * g_.Gk(p<0>(w) + p<2>(w) - p<1>(w), back_to_bz( p<0>(k)-p<0>(qpr) ), back_to_bz(p<1>(k)-p<1>(qpr) )) );//std::pow(k2mesh.size(),-2.0) ;
           }
           for(Mesh<FERMI>::iterator f(fmesh); !f.is_end(); f.increment() )
               std::cout<< f.val() << " " << p<0>(k) << " " << p<1>(k)<<" " << s_.Gk(f.r_index() , p<0>(k) , p<1>(k))[0][0] << " " << s_.Gk(f.r_index(), p<0>(k) , p<1>(k))[0][1] << "\n";
  
        }
  
      } else {
        double dos4pm2 = std::pow(ek.dos[4], -2.0);
      for ( mesh_product< Lmesh, Lmesh>::iterator k(k2mesh); !k.is_end(); k.increment() ) {
              for ( mesh_product< Lmesh, Lmesh>::iterator q(q2mesh); !q.is_end(); q.increment() )
                  for ( mesh_product< Lmesh, Lmesh>::iterator qpr(qprim2mesh); !qpr.is_end(); qpr.increment() ){
                     for (mesh_product<Mesh<FERMI> , Mesh<FERMI>, Mesh<FERMI> >::iterator w(iw_mesh); !w.is_end(); w.increment() ) {
                      for (int i=0; i<2; ++i)
                          for ( int j =0 ; j<2; ++j){
                              A[i][j] = g_.Gk( p<1>(w), back_to_bz(p<0>(k)-p<0>(q)) ,back_to_bz( p<1>(k)-p<1>(q)) )[i][j]; 
                              B[i][j] = g_.Gk(p<2>(w), back_to_bz(p<0>(k)-p<0>(q)-p<0>(qpr)), back_to_bz(p<1>(k)-p<1>(q)-p<1>(qpr) ) )[i][j] ; 
                              C[i][j] = g_.Gk( p<0>(w)+p<2>(w)-p<1>(w), back_to_bz(p<0>(k)-p<0>(qpr)) ,back_to_bz( p<1>(k)-p<1>(qpr)) )[i][j];
                          }
                       temp3 =  A[0][0] * B[0][0] - B[1][0] * A[0][1];
                       temp4 = -B[0][1] * A[0][0] + B[1][1] * A[0][1];
                       kern = ( Ker( abs(p<0>(w)-p<1>(w) ) ) * Ker( abs( p<1>(w) - p<2>(w) ) ) ) * dos4pm2;
                       s_.Gk(p<0>(w), p<0>(k) , p<1>(k))[0][0]  += (C[0][0] * temp3 + C[1][0] * temp4 ) * kern * std::pow(k2mesh.size(),-2.0);
                       s_.Gk(p<0>(w), p<0>(k) , p<1>(k))[0][1]  += (C[0][1] * temp3 + C[1][1] * temp4 ) * kern * std::pow(k2mesh.size(),-2.0);
                      }
//                      std::cout << " ================== \n";
//                      for(Mesh<FERMI>::iterator f(fmesh); !f.is_end(); f.increment() )
//                        
//                          std::cout<< f.val() << " k= " << p<0>(k) << " " << p<1>(k)<<" "<< " q= " << p<0>(q) << " " << p<1>(q)<<" "" qpr= " << p<0>(qpr) << " " << p<1>(qpr)<<" " << s_.Gk(f.r_index() , p<0>(k) , p<1>(k))[0][0] << " " << s_.Gk(f.r_index(), p<0>(k) , p<1>(k))[0][1] << "\n";
                }
//             for(Mesh<FERMI>::iterator f(fmesh); !f.is_end(); f.increment() )
//                          std::cout<< f.val() << " k= " << p<0>(k) << " " << p<1>(k)<<" " << s_.Gk(f.r_index() , p<0>(k) , p<1>(k))[0][0] << " " << s_.Gk(f.r_index(), p<0>(k) , p<1>(k))[0][1] << "\n";
          }
        }
       
        t  = clock() - t;
        printf ("Frequency summation for second diagram (%f seconds).\n",((float)t)/CLOCKS_PER_SEC);

        Delta delta;
        double sigma_ = 0.01;
        double Tp2 = pow(T,2.0);
//std::cout << " llllll\n";
        for (Mesh<FERMI>::iterator f(fmesh); !f.is_end(); f.increment() ){
           for ( mesh_product< Lmesh, Lmesh>::iterator k(k2mesh); !k.is_end(); k.increment() ) 
               S2( f.r_index() ) += delta( (ek( k )[4]-ek.Ef)/sigma_ )/sigma_/k2mesh.size()/ek.dos[4] * s_.Gk( f.r_index(), p<0>(k), p<1>(k) );
           S2( f.r_index() )[1][1] = -std::conj( S2( f.r_index() )[0][0] ); S2( f.r_index() )[1][0] = S2( f.r_index() )[0][1];
           S( f.r_index() ) += s3 * ( Tp2 * S2( f.r_index() ) * s3);
        }
    }
     
//std:: cout << " ffff2 \n" ;
return;
}
        
  
void iterate(){
     std::complex<double> alpha( 0.05,0);
     GF< Mesh<FERMI>, std::complex<double>, multi  > S_old( S );
     Mesh<FERMI> fmesh3(T, 3*Niw );
     std::cout<<"*  -----------------------------------------------------------  *\n";
     std::cout<<"Superconducting temperature = "<<T<<", Nw = "<<Niw<<", muc = "<<muc<<", mucn = "<<mucn<<"\n";
     std::cout<<"\n";std::cout<<"\n";
     std::ofstream of; 
     double phi_old(0.0),error(1.0),zi_old(0);
     for( iter=0; iter<maxiter; ++iter){

        // mixing , simple, broyden what ever...
        for (Mesh<FERMI>::iterator f(fmesh); !f.is_end(); f.increment() ){
            S(f.r_index()) = alpha*S_old(f.r_index()) + (1.0-alpha)*S(f.r_index());
        }
     

           
        S_old = S;
                std::cout<<"*  -----------------------------------------------------------  *\n";
        std::cout<<"*  --- Iteration          "<<iter<<"                     -----  *\n";
        std::cout << " T = " << T/8.6173427909e-05 << "  phcut = " << phcut <<"\n";
    //  eval_Sigma_direct();
        eval_Sigma();
        phi_old = 0.5*( S_old(fmesh3(fmesh3.size()/2))[0][1] + S_old(fmesh3(fmesh3.size()/2))[1][0] ).real()/(1.0 - S_old(fmesh3(fmesh3.size()/2))[0][0].imag()/fmesh3(fmesh3.size()/2));
        zi_old = 1.0 - S_old(fmesh3(fmesh3.size()/2))[0][0].imag()/fmesh3(fmesh3.size()/2);
        //
        of.open(std::string("Delta_")+basename+std::to_string(iter));
        for (Mesh<FERMI>::iterator f(fmesh3); !f.is_end(); f.increment() ){
            zi(f) = 1.0 - S(f)[0][0].imag()/f.val();         
            phi(f) = 0.5*( S(f)[0][1] + S(f)[1][0] ).real()*std::pow(zi(f),-1.0); 
            chi(f)  = S(f)[0][0].real()*std::pow(zi(f),-1.0);
            of<< f.val() << " " << phi(f) << " " <<  zi(f)<<" " << chi(f) << "\n"; 
            
        }
        error = std::abs((phi(fmesh3.size()/2)-phi_old)/phi_old); 
        of.close();
        std::cout << "\n";
        std::cout<< "       w0       " << "        delta[0](phi/zi)         " << "     zi[0]          error(%)\n "; 
        std::cout<<"    "<<fmesh3(fmesh3.size()/2)<<"               " << phi(fmesh3.size()/2) << "               "<< zi(fmesh3.size()/2) << "          "<<error << "\n";
        std::cout<<"\n";std::cout<<"\n";
        std::cout<<"*  -----------------------------------------------------------  *\n";
        if ( error < 10e-5 )
           break;
     }
//std::cout << phi ;
//std::cout << " \n ";
//std::cout << zi ;
return;
}

private:
int maxiter;
GF< Mesh<FERMI>, double, scalar > phi, chi, zi;
Kernel<Mesh<BOSON>, double > Ker;
matrix< std::complex<double> > s0, s1, s2, s3;
GF<Mesh<FERMI>, std::complex<double>, multi > S;
Mesh<FERMI> fmesh;
Dmesh rmesh;
Lmesh kmesh;
bool second_diag;
dispersion<DISPSRC>  ek;
int Niw,iter,Nr;
double T,wscut,muc,mucn,Ef,phcut,w2,wn;
std::string basename;
};

#include <stdlib.h>
//
int main(int argc, char **argv)
{

if (atof(*(argv+1)) != 0 ){
double sctemp_min= atof(*(argv+1));// Kelvinsctempsctemp
double sctemp_max=  atof(*(argv+2));
double Ef = atof(*(argv+3));
//int N_tempr = 10;
double sctemp(0.0);
bool vertex_flag(false);
double muc= atof(*(argv+4)); 
double wscut=atof(*(argv+5)); // ev        
std::string basename;
double phcut =  atof(*(argv+7));
std::string sn(*(argv+8));
double N_tempr = atof(*(argv+9));
double sigma = atof(*(argv+10));
int Nk = atoi(*(argv+11));
int Nr = atoi(*(argv+12));
eliashberg<second_diagram, fromfile >* e1; 

for ( int a=0; a< N_tempr; ++a) {
  sctemp =  sctemp_min + a*(sctemp_max-sctemp_min)/N_tempr;
  if ( atof(*(argv+6)) == 1 ){
   vertex_flag = true;
   basename =  std::string("mucNr")+std::to_string(Nr)+ std::string("Nk") + std::to_string(Nk)+std::string("sigma") +std::to_string( sigma ) + std::string(sn.data()) + std::string("_newmucker") + std::string("T") +  std::to_string(sctemp) + std::string("Ef") + std::string( *(argv+3) ) + std::string("muc") + std::string( *(argv+4) )+std::string("wscut") + std::string(*(argv+5))+ std::string("phcut")+ std::to_string(phcut) + std::string("1ndvertex");
} else {
  basename =  std::string("mucNr")+std::to_string(Nr)+ std::string("Nk") + std::to_string(Nk)+std::string("sigma") +std::to_string( sigma ) + std::string(sn.data()) + std::string("_newmucker") + std::string("T") +  std::to_string( sctemp ) + std::string("Ef") + std::string( *(argv+3) ) + std::string("muc") + std::string( *(argv+4) )+std::string("wscut") + std::string(*(argv+5)) + std::string("phcut")+ std::to_string(phcut) + std::string("novertex");
}
  sctemp*=kelvin2eV;
  int Nw = wscut/(pi*sctemp) -1. ;
  if (Nw - 2*(Nw/2) == 1) Nw--;
  double init_gap(0);
  string s ("/home/davoud-i7/eliashberg/band.eig_");
  s = s + std::to_string(Nk) + "_3";
  std::cout << "Read dispersion from : " << s << "\n";
  if ( sn == "sup")
     init_gap = 0.005;
  std::cout << " The initial state is : " <<sn.data() <<" with initial gap: "<<init_gap<<"\n"; 
  eliashberg<second_diagram, fromfile > *e2= new eliashberg<second_diagram, fromfile >(sctemp, muc, Nw,Nk, init_gap,
                                        vertex_flag , std::string("pb11_nodipvac20ecut90k2020_tr1.0d-15_q16.a2f.01"),
                                        dispersion<fromfile>(s, 8, Ef, sigma, Nk),phcut, basename, Nr);

e2->iterate();
delete e2; 
}

} else {

std::cout << " el.x tmin tmax Ef muc wscut vertex phcut sup/nor temps_steps sigma Nk Nr\n";
}
return 0;



}

