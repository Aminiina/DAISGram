#include <iostream>
#include <string>
#include <random>
#include <math.h>
#include <fstream>

#include "dais_exc.h"
#include "tensor.h"

#define PI 3.141592654
#define FLT_MAX 3.402823466e+38F /* max value */
#define FLT_MIN 1.175494351e-38F /* min positive value */

using namespace std;

void Tensor::init_random(float mean, float std){
    if(data){

        std::default_random_engine generator;
        std::normal_distribution<float> distribution(mean,std);

        for(int i=0;i<r;i++){
            for(int j=0;j<c;j++){
                for(int k=0;k<d;k++){
                    this->operator()(i,j,k)= distribution(generator);
                }
            }
        }

    }else{
        throw(tensor_not_initialized());
    }
}

Tensor::Tensor(){

    data = nullptr;
    r = 0;
    c = 0;
    d = 0;

}
Tensor::~Tensor(){

    for(int i=0;i<r;i++){

        for(int j=0;j<c;j++){
            delete[] data[i][j];
        }
        delete[] data[i];
    }
    delete[] data;
}

Tensor::Tensor(int r, int c, int d, float v){

    if(r<=0 || c<=0 || d<=0) throw index_out_of_bound();

    this->r = r;
    this->c = c;
    this->d = d;

    data = new float**[r];

    for(int i=0;i<r;i++){
        data[i] = new float*[c];
        for(int j=0;j<c;j++){
            data[i][j] = new float[d];

        }
    }

    for(int i=0;i<r;i++){
        for(int j=0;j<c;j++){
            for(int k=0;k<d;k++){
                data[i][j][k] = v;

            }
        }
    }

}



float Tensor::operator()(int i, int j, int k) const{

    if (i < 0 || j < 0 || k < 0 || i >= r || j >= c || k >= d ){
        throw(index_out_of_bound());
    }
    return data[i][j][k];
}

float& Tensor::operator()(int i, int j, int k){

    if (i < 0 || j < 0 || k < 0 || i >= r || j >= c || k >= d ){
        throw index_out_of_bound();
    }
    return data[i][j][k]; // data(i,j,k)

}

 Tensor::Tensor(const Tensor& that){

    if(!that.data) throw tensor_not_initialized();

    this->r = that.r;
    this->c = that.c;
    this->d = that.d;

    this->data = new float**[r];

    for(int i=0;i<this->r;i++){
        this->data[i] = new float*[this->c];
        for(int j=0;j<this->c;j++){
            this->data[i][j] = new float[this->d];

        }
    }

    for(int i=0;i<this->r;i++){
        for(int j=0;j<this->c;j++){
            for(int k=0;k<this->d;k++){
                this->data[i][j][k] = that.data[i][j][k];
            }
        }

    }
}

bool Tensor::operator==(const Tensor& rhs) const{
    if(!data || !rhs.data) throw tensor_not_initialized();
    if(this->r != rhs.r || this->c != rhs.c || this->d != rhs.d){
        throw dimension_mismatch();
    }else {
        bool risp=false;

            for(int i=0;i<r;i++){
                for(int j=0;j<c;j++){
                    for(int k=0;k<d;k++){
                        if(abs(this->data[i][j][k]-rhs.data[i][j][k])<EPSILON){
                            risp=true;
                        }
                    }
                }
            }
        return risp;
    }
}

Tensor Tensor::operator-(const Tensor &rhs)const{
    if(!data || !rhs.data) throw tensor_not_initialized();
    if(this->r != rhs.r || this->c != rhs.c || this->d != rhs.d){
        throw dimension_mismatch();
    }
    else{
        Tensor newTensor{r,c,d,0.0};

        for(int i=0;i<r;i++){
            for(int j=0;j<c;j++){
                for(int k=0;k<d;k++){
                    newTensor(i,j,k)= data[i][j][k]-rhs(i,j,k);
                }
            }
        }
        return newTensor;
    }
}

Tensor Tensor::operator +(const Tensor &rhs)const{
    if(!data || !rhs.data) throw tensor_not_initialized();
    if(this->r != rhs.r || this->c != rhs.c || this->d != rhs.d){
        throw dimension_mismatch();
    }
    else{
        Tensor newTensor=Tensor(r,c,d,0.0);
        for(int i=0;i<r;i++){
            for(int j=0;j<c;j++){
                for(int k=0;k<d;k++){
                    newTensor.data[i][j][k]=this->data[i][j][k]+rhs.data[i][j][k];
                }
            }
        }
        return newTensor;
    }
}

Tensor Tensor::operator*(const Tensor &rhs)const{
    if(!data || !rhs.data) throw tensor_not_initialized();
    if(this->r != rhs.r || this->c != rhs.c || this->d != rhs.d){
        throw dimension_mismatch();
    }
    else{
        Tensor newTensor=Tensor(r,c,d,0.0);
        for(int i=0;i<r;i++){
            for(int j=0;j<c;j++){
                for(int k=0;k<d;k++){
                    newTensor.data[i][j][k]=this->data[i][j][k]*rhs.data[i][j][k];
                }
            }
        }
        return newTensor;
    }
}

Tensor Tensor::operator/(const Tensor &rhs)const{
    if(!data || !rhs.data) throw tensor_not_initialized();
    if(this->r != rhs.r || this->c != rhs.c || this->d != rhs.d){
        throw dimension_mismatch();
    }
    else{
        Tensor newTensor=Tensor(r,c,d,0.0);
        for(int i=0;i<r;i++){
            for(int j=0;j<c;j++){
                for(int k=0;k<d;k++){
                    newTensor.data[i][j][k]=this->data[i][j][k]/rhs.data[i][j][k];
                }
            }
        }
        return newTensor;
    }
}

Tensor Tensor::operator-(const float &rhs)const {
    if (!data) throw tensor_not_initialized();
    else {
        Tensor newTensor = Tensor(r, c, d, 0.0);
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                for (int k = 0; k < d; k++) {
                    newTensor.data[i][j][k] = this->data[i][j][k] - rhs;
                }
            }
        }
        return newTensor;
    }
}

Tensor Tensor::operator+(const float &rhs)const{
    if(!data) throw tensor_not_initialized();
    else {
        Tensor newTensor = Tensor(r, c, d, 0.0);
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                for (int k = 0; k < d; k++) {
                    newTensor.data[i][j][k] = this->data[i][j][k] + rhs;
                }
            }
        }
        return newTensor;
    }
}

Tensor Tensor::operator*(const float &rhs)const{
    if(!data) throw tensor_not_initialized();
    else {
        Tensor newTensor = Tensor(r, c, d, 0.0);
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                for (int k = 0; k < d; k++) {
                    newTensor.data[i][j][k] = this->data[i][j][k] * rhs;
                }
            }
        }
        return newTensor;
    }
}

Tensor Tensor::operator/(const float &rhs)const{
    if(!data) throw tensor_not_initialized();
    else {
        Tensor newTensor = Tensor(r, c, d, 0.0);
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                for (int k = 0; k < d; k++) {
                    newTensor.data[i][j][k] = this->data[i][j][k] / rhs;
                }
            }
        }
        return newTensor;
    }
}

Tensor& Tensor:: operator=(const Tensor &other){

    if (data){
        for(int i=0;i<r;i++){
            for(int j=0;j<c;j++){
                delete[] data[i][j];
            }
            delete[] data[i];
        }
        delete[] data;
    }

    this->r = other.r;
    this->c = other.c;
    this->d = other.d;

    this->data = new float**[r];

    for(int i=0;i<this->r;i++){
        this->data[i] = new float*[this->c];
        for(int j=0;j<this->c;j++){
            this->data[i][j] = new float[this->d];

        }
    }

    for(int i=0;i<this->r;i++){
        for(int j=0;j<this->c;j++){
            for(int k=0;k<this->d;k++){
                this->data[i][j][k] = other.data[i][j][k];
            }
        }

    }
    return *this;
}

void Tensor::init(int r, int c, int d, float v){

    if(r < 0 || c < 0 || d < 0) throw index_out_of_bound();

    if (data){
        for(int i=0;i<r;i++){
            for(int j=0;j<c;j++){
                delete[] data[i][j];
            }
            delete[] data[i];
        }
        delete[] data;
    }

    *this = Tensor(r,c,d,v);

}

void Tensor::clamp(float low, float high){

    if(data){
        for (int i=0; i<r; i++){
            for (int j=0;j<c;j++){
                for (int k=0;k<d;k++){
                    if (this->operator()(i,j,k)<low){
                        this->operator()(i,j,k)=low;
                    }
                    else if(this->operator()(i,j,k)>high){
                        this->operator()(i,j,k)=high;
                    }
                }
            }
        }
    } else throw tensor_not_initialized();
}

void Tensor::rescale(float new_max){

    if (data){
        for (int k=0;k<d;k++){
            float min = getMin(k);
            float max = getMax(k);
            for (int i=0; i<r; i++){
                for (int j=0;j<c;j++){
                    if (max - min) (*this)(i,j,k)=(((*this)(i,j,k) - min)/(max-min))*new_max;
                    else (*this)(i,j,k)= 0;
                }
            }
        }
    } else throw tensor_not_initialized();
}

Tensor Tensor::padding(int pad_h, int pad_w)const{
    if (!data) throw tensor_not_initialized();

    int new_r=r+2*pad_h;
    int new_c=c+2*pad_w;

    Tensor newTensor=Tensor(new_r,new_c,d,0.0);

    for(int i=0;i<this->r;i++){
        for(int j=0;j<this->c;j++){
            for(int k=0;k<this->d;k++){
                newTensor.data[i+pad_h][j+pad_w][k] = this->data[i][j][k];
            }
        }
    }
    return newTensor;
}

Tensor Tensor::subset(unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end, unsigned int depth_start, unsigned int depth_end)const{

    if (!data) throw tensor_not_initialized();
    if (row_end >= r || col_end >= c || depth_end >= d || row_end < 0 || col_end < 0 || depth_end < 0 || row_start >= r || col_start >= c || depth_start >= d || row_start < 0 || col_start < 0 || depth_start < 0 ) throw(index_out_of_bound());

    int i,j,k=0;
    i=row_end-row_start;
    j=col_end-col_start;
    k=depth_end-depth_start;

    Tensor newTensor=Tensor(i,j,k,0.0);
    for (int a=0; a<i; a++ ){
        for(int b=0; b<j; b++){
            for(int c=0;c<k;c++){
                newTensor.data[a][b][c]=this->data[a+row_start][b+col_start][c+depth_start];
            }
        }
    }
    return newTensor;
}

Tensor Tensor::concat(const Tensor &rhs, int axis)const{
    if (!data || !rhs.data) throw tensor_not_initialized();
    if (axis < 0 || axis > 2) throw unknown_exception();

    if (axis==0 && this->c==rhs.c && this->d==rhs.d){

        int new_r = this->r + rhs.r;
        Tensor newTensor (new_r,c,d,0.0);
        for (int i = 0; i < new_r; i++){
            for (int j = 0; j < newTensor.c; j++){
                for (int k = 0; k < newTensor.d; k++){
                    if (i < this->r ) newTensor.data[i][j][k] =this->data[i][j][k];
                    else newTensor.data[i][j][k] =rhs.data[i-this->r][j][k];
                }
            }
        }
        return newTensor;
    }
    if (axis==1 && this->r==rhs.r && this->d==rhs.d){
        int new_c = this->c + rhs.c;
        Tensor newTensor (r,new_c,d,0.0);
        for (int i = 0; i < newTensor.r; i++){
            for (int j = 0; j < new_c; j++){
                for (int k = 0; k < newTensor.d; k++){
                    if (j < this->c ) newTensor.data[i][j][k] =this->data[i][j][k];
                    else newTensor.data[i][j][k] =rhs.data[i][j-this->c][k];
                }
            }
        }
        return newTensor;
    }
    if (axis==2 && this->r==rhs.r && this->c==rhs.c){
        int new_d = this->d + rhs.d;
        Tensor newTensor (r,c,new_d,0.0);
        for (int i = 0; i < newTensor.r; i++){
            for (int j = 0; j < newTensor.c; j++){
                for (int k = 0; k < new_d; k++){
                    if (k < this->d ) newTensor.data[i][j][k] =this->data[i][j][k];
                    else newTensor.data[i][j][k] =rhs.data[i][j][k-this->d];
                }
            }
        }
        return newTensor;
    }
    else{
        throw(concat_wrong_dimension());
    }
}

Tensor Tensor::convolve(const Tensor &f)const{

    if (!data || !f.data) throw tensor_not_initialized();

    if(this->depth()!=f.depth()){
        throw(dimension_mismatch());
    }
    if(f.rows()%2==0 or f.cols()%2==0){
        throw(filter_odd_dimensions());
    }

    Tensor newTensor(r,c,d);
    Tensor padded = this->padding((f.r-1)/2,(f.c-1)/2);
    float temp = 0;

    for (int k = 0; k < padded.d; k++){
        for (int i = 0; i < padded.r-f.r+1; i++){
            for (int j = 0; j < padded.c-f.c+1; j++){
                    for (int b = 0; b < f.c; b++){
                        for (int a = 0; a < f.r; a++){
                            temp += padded.data[i+a][j+b][k]*f.data[a][b][k];
                        }
                    }
                newTensor.data[i][j][k] = temp;
                temp = 0;
            }
        }
    }
    return newTensor;
}

int Tensor::rows()const{
    if (!data) throw tensor_not_initialized();
    else return r;
}

int Tensor::cols()const{
    if (!data) throw tensor_not_initialized();
    else return c;
}

int Tensor::depth()const{
    if (!data) throw tensor_not_initialized();
    else return d;
}

float Tensor::getMin(int k) const{

    if (!data) throw tensor_not_initialized();

    if (k >= d or k < 0) throw index_out_of_bound();

    float min{(*this)(0,0,k)};

    for (int i=0; i<r; i++){
        for (int j=0; j<c; j++){
            if (this->data[i][j][k] < min) min = this->data[i][j][k];
        }
    }
    return min;
}

float Tensor::getMax(int k)const{

    if (!data) throw tensor_not_initialized();

    if (k >= d or k < 0) throw index_out_of_bound();

    float max = this->data[0][0][k];
    for (int i=0; i<r; i++){
        for (int j=0; j<c; j++){
            if (this->data[i][j][k] > max) max = this->data[i][j][k];

        }
    }
    return max;
}

void Tensor::showSize()const{
    if(data){
        std::cout<<"rows"<<this->rows()<<"columns"<<this->cols()<<"depth"<<this->depth();
    } else throw tensor_not_initialized();
}

ostream& operator<< (ostream& stream, const Tensor & obj){
    if(!obj.data) throw tensor_not_initialized();

    for(int k=0;k<obj.d;k++){
        for(int i=0;i<obj.r;i++){
            stream << "( ";
            for(int j=0;j<obj.c;j++){
                stream << obj(i,j,k) << " ";
            }
            stream << ")";
            stream << "\n";
        }
        stream << "\n";
    }
    return stream;
}

void Tensor::read_file(string filename){

    if (filename.empty()) throw unable_to_read_file();

    ifstream is{filename};

    is >> r;
    is >> c;
    is >> d;

    init(r,c,d);

    if(data){
        for(int k=0;k<d;k++){
            for(int i=0;i<r;i++){
                for(int j=0;j<c;j++){
                    is >> this->operator()(i,j,k);
                }
            }
        }
    }
}

void Tensor::write_file(string filename){

    ofstream os{filename};

    os << r << "\n";
    os << c << "\n";
    os << d << "\n";

    if(data){
        for(int k=0;k<d;k++){
            for(int i=0;i<r;i++){
                for(int j=0;j<c;j++){
                    os << this->operator()(i,j,k) << "\n";
                }
            }
        }
    }
}

