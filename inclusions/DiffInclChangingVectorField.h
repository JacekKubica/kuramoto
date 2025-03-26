//#include "DiffInclVectorField.h"

//class DiffInclChangingVectorField : public DiffInclVectorField {
//public:
    //// TODO from base class 
    //typedef capd::interval ScalarType;
    //typedef capd::IVector VectorType;
    //typedef capd::IMatrix MatrixType;
    //typedef capd::IMap MapType;
    //DiffInclChangingVectorField(MapType &selector, MapType &perturbation, MapType &full_system,
                                //MapType &cut_out_system, size_t m) 
        //: DiffInclVectorField(selector, perturbation), full_system(full_system),
          //cut_out_system(cut_out_system), m(m) {}
    //VectorType operator()(const ScalarType t, const VectorType &x) {
        //return full_system(t, x);
    //}
    //void setSelector(const VectorType &W2) override {
        //auto perturbs = cut_out_system(W2);
        //for(size_t i = 0; i < m; ++i) {
            //getPerturbation().setParameter(i, perturbs[i]);
            //getSelector().setParameter(i, -perturbs[i]);
        //}
    //}
    //VectorType cutDimension(const VectorType W2) {
        //VectorType W2_(m);
        //for(size_t i = 0; i < m; ++i) {
            //W2_[i] = W2[i];
        //}
        //return W2_;
    //}
//private:
    //MapType full_system;
    //MapType cut_out_system;
    //size_t m;
//};

