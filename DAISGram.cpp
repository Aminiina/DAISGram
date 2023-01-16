#include <iostream>
#include <string>

#include "dais_exc.h"
#include "tensor.h"
#include "libbmp.h"
#include "DAISGram.h"

using namespace std;

/**
 * Load a bitmap from file
 *
 * @param filename String containing the path of the file
 */
void DAISGram::load_image(string filename){
    BmpImg img = BmpImg();

    img.read(filename.c_str());

    const int h = img.get_height();
    const int w = img.get_width();

    data = Tensor(h, w, 3, 0.0);

    for(int i=0;i<img.get_height();i++){
        for(int j=0;j<img.get_width();j++){ 
            data(i,j,0) = (float) img.red_at(j,i);
            data(i,j,1) = (float) img.green_at(j,i);    
            data(i,j,2) = (float) img.blue_at(j,i);   
        }                
    }
}


/**
 * Save a DAISGram object to a bitmap file.
 * 
 * Data is clamped to 0,255 before saving it.
 *
 * @param filename String containing the path where to store the image.
 */
void DAISGram::save_image(string filename){

    data.clamp(0,255);

    BmpImg img = BmpImg(getCols(), getRows());

    img.init(getCols(), getRows());

    for(int i=0;i<getRows();i++){
        for(int j=0;j<getCols();j++){
            img.set_pixel(j,i,(unsigned char) data(i,j,0),(unsigned char) data(i,j,1),(unsigned char) data(i,j,2));                   
        }                
    }

    img.write(filename);

}

DAISGram::DAISGram(){}

DAISGram::~DAISGram(){}

int DAISGram::getRows(){
    return data.rows();
}

int DAISGram::getCols(){
    return data.cols();
}

int DAISGram::getDepth(){
    return data.depth();
}

DAISGram DAISGram::brighten(float bright){

    if(bright<0) throw unknown_exception();
    if(getRows() < 0 || getCols() < 0 || getDepth() < 0) throw index_out_of_bound();

    DAISGram newdais(*this);
    newdais.data=this->data;
    for(int i=0;i<getRows();i++){
        for(int j=0;j<this->getCols();j++){
            for(int k=0;k<this->getDepth();k++){
                newdais.data.operator()(i,j,k)+=bright;
            }
        }
    }
    newdais.data.clamp(0,255);
    return newdais;

}

DAISGram DAISGram::grayscale(){

    if(getRows() < 0 || getCols() < 0 || getDepth() < 0) throw index_out_of_bound();
    if (getDepth()>3) throw unknown_exception();

    DAISGram newdais(*this);
    int tot=0;
    float somma=0;
    float media = 0;
    for(int i = 0; i < getRows(); i++) {
        for (int j = 0; j < getCols(); j++) {
            for (int k = 0; k < getDepth(); k++) {
                somma += data.operator()(i, j, k);
                tot ++;
            }
            media = somma / tot;
            for (int k = 0; k < getDepth(); k++) {
                newdais.data(i, j, k) = media;
            }
            somma = 0;
            tot = 0;
        }
    }
    return newdais;
}

DAISGram DAISGram::warhol(){
    if(getRows() < 0 || getCols() < 0 || getDepth() < 0) throw index_out_of_bound();
    if (getDepth()>3) throw unknown_exception();

    Tensor newTensor(getRows()*2,getCols()*2,getDepth());
    DAISGram newdais;
    newdais.data = newTensor;

    for (int i = 0; i < this->getRows(); i++){//top left
        for (int j = 0; j < this->getCols(); j++){
            for (int k = 0; k < this->getDepth(); k++){
                newdais.data(i,j,k) = this->data(i,j,k);
            }
        }
    }

    for (int k = 0; k < this->getDepth(); k++){//top right
        for (int i = 0; i < this->getRows(); i++){
            for (int j = this->getCols(); j < newdais.getCols(); j++){
                if (k == 0) newdais.data(i,j,0) = this->data(i,j-this->getCols(),1);
                if (k == 1) newdais.data(i,j,1) = this->data(i,j-this->getCols(),0);
                if (k == 2) newdais.data(i,j,k) = this->data(i,j-this->getCols(),k);
            }
        }
    }

    for (int i = this->getRows(); i < newdais.getRows(); i++){//bottom left
        for (int j = 0; j < this->getCols(); j++){
            for (int k = 0; k < this->getDepth(); k++){
                if (k == 2) newdais.data(i,j,k) = this->data(i - this->getRows(),j,1);
                if (k == 1) newdais.data(i,j,k) = this->data(i - this->getRows(),j,2);
                if (k == 0) newdais.data(i,j,k) = this->data(i - this->getRows(),j,k);
            }
        }
    }

    for (int i = this->getRows(); i < newdais.getRows(); i++){//bottom right
        for (int j = this->getCols(); j < newdais.getCols(); j++){
            for (int k = 0; k < this->getDepth(); k++){
                if (k == 2) newdais.data(i,j,k) = this->data(i - this->getRows(),j-this->getCols(),0);
                if (k == 0) newdais.data(i,j,k) = this->data(i - this->getRows(),j-this->getCols(),2);
                if (k == 1)  newdais.data(i,j,k) = this->data(i - this->getRows(),j-this->getCols(),k);
            }
        }
    }

    return newdais;
}

DAISGram DAISGram::sharpen() {
    if(getRows() < 0 || getCols() < 0 || getDepth() < 0) throw index_out_of_bound();
    if (getDepth()>3) throw unknown_exception();

    Tensor filter = Tensor(3, 3, getDepth());

    for (int k = 0; k < getDepth(); k++){
        filter(0, 0, k) = 0;
        filter(0, 1, k) = -1;
        filter(0, 2, k) = 0;
        filter(1, 0, k) = -1;
        filter(1, 1, k) = 5;
        filter(1, 2, k) = -1;
        filter(2, 0, k) = 0;
        filter(2, 1, k) = -1;
        filter(2, 2, k) = 0;
    }

    DAISGram newdais(*this);
    newdais.data=data.convolve(filter);
    newdais.data.clamp(0,255);
    return newdais;
}

DAISGram DAISGram::emboss(){
    if(getRows() < 0 || getCols() < 0 || getDepth() < 0) throw index_out_of_bound();
    if (getDepth()>3) throw unknown_exception();

    Tensor filter=Tensor(3,3,getDepth());
    for (int k = 0; k < getDepth(); k++) {
        filter(0, 0, k) = -2;
        filter(0, 1, k) = -1;
        filter(0, 2, k) = 0;
        filter(1, 0, k) = -1;
        filter(1, 1, k) = 1;
        filter(1, 2, k) = 1;
        filter(2, 0, k) = 0;
        filter(2, 1, k) = 1;
        filter(2, 2, k) = 2;
    }

    DAISGram newdais(*this);
    newdais.data=data.convolve(filter);
    newdais.data.clamp(0,255);
    return newdais;

}

DAISGram DAISGram::smooth(int h){
    if(getRows() < 0 || getCols() < 0 || getDepth() < 0) throw index_out_of_bound();
    if (getDepth()>3) throw unknown_exception();
    if(h<0) throw unknown_exception();

    float c=1.0/((float)(h*h));
    Tensor filter(h,h,getDepth(),c);
    DAISGram newdais(*this);
    newdais.data=data.convolve(filter);
    return newdais;
}

DAISGram DAISGram::edge(){
    if(getRows() < 0 || getCols() < 0 || getDepth() < 0) throw index_out_of_bound();
    if (getDepth()>3) throw unknown_exception();

    Tensor filter=Tensor(3,3,getDepth());
    for (int k = 0; k < getDepth(); k++) {
        filter(0, 0, k) = -1;
        filter(0, 1, k) = -1;
        filter(0, 2, k) = -1;
        filter(1, 0, k) = -1;
        filter(1, 1, k) = 8;
        filter(1, 2, k) = -1;
        filter(2, 0, k) = -1;
        filter(2, 1, k) = -1;
        filter(2, 2, k) = -1;
    }

    DAISGram newdais(*this);
    DAISGram image_gray=this->grayscale();
    newdais.data=image_gray.data.convolve(filter);
    newdais.data.clamp(0,255);
    return newdais;

}

DAISGram DAISGram::blend(const DAISGram & rhs, float alpha) {
    //if(getRows() < 0 || getCols() < 0 || getDepth() < 0) throw index_out_of_bound();
    if (alpha<0 || alpha > 1) throw unknown_exception();

    DAISGram newdais(*this);
    newdais.data = data*alpha + rhs.data*(1-alpha);
    return newdais;
}

DAISGram DAISGram::greenscreen(DAISGram & bkg, int rgb[], float threshold[]){

    if (getRows() != bkg.getRows() || getCols() != bkg.getCols() || getDepth()!= bkg.getDepth()) throw dimension_mismatch();

    DAISGram newdais(*this);
    int count = 0;

    for (int i = 0; i < this->getRows(); i++){
        for (int j = 0; j < this->getCols(); j++){
            for (int k = 0; k < this->getDepth(); k++){
                if (newdais.data(i,j,k) >= (rgb[k]-threshold[k]) && newdais.data(i,j,k) <= (rgb[k]+threshold[k])) count++;
            }
            if (count == 3) {
                for (int k = 0; k < this->getDepth(); k++) {
                    newdais.data(i, j, k) = bkg.data(i, j, k);
                }
            }
            count = 0;
        }
    }
    return newdais;
}

DAISGram DAISGram::equalize(){
    throw method_not_implemented();

}

/**
 * Generate Random Image
 * 
 * Generate a random image from nois
 * 
 * @param h height of the image
 * @param w width of the image
 * @param d number of channels
 * @return returns a new DAISGram containing the generated image.
 */  
void DAISGram::generate_random(int h, int w, int d){
    data = Tensor(h,w,d,0.0);
    data.init_random(128,50);
    data.rescale(255);
}