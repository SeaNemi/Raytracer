#include<iostream>
#include <eigen3/Eigen/Dense>
#include <string>
#include <getopt.h>
#include "Raytracer.h"

int main(int argc, char* argv[]){
    std::cout << "Ray Tracer!" << std::endl;
    int c;
    double aperture = 0.0;
    int samples = 1;
    int maxraydepth = 5;
    bool color = false;
    bool phong = false;
    bool stratified = false;
    bool reflections = true;
    bool shadows = true;
    bool verbose = false;
    int mode = 3;
    while ((c = getopt(argc, argv, "a:s:d:cpjm:v")) != -1) {
        switch(c) {
        case 'a':
        aperture = atof(optarg);
        break;
        case 's':
        samples = atoi(optarg);
        break;
        case 'j':
        stratified = true;
        break;
        case 'c':
        color = true;
        break;
        case 'p':
        phong = true;
        break;
        case 'm':
        mode = atoi(optarg);
        break;
        case 'd':
        maxraydepth = atoi(optarg);
        break;
        case 'v':
        verbose = true;
        break;
        default:
        abort();
        }
    }

    if (argc-optind != 2) {
        std::cout<<"usage: trace input.nff output.ppm"<<std::endl;
        for (int i=0; i<argc; i++) std::cout<<argv[i]<<std::endl;
        exit(0);
    }	

    switch(mode) {
        case 0:
            reflections = shadows = false;
            break;
        case 1:
            reflections = false;
            break;
        case 2:
            shadows = false;
            break;
    }

    Scene newScene(argv[optind++]);
    if(newScene.parse()){
        //scene then sent to the ray tracer
        Raytracer tracer(newScene);
        tracer.aperture = aperture;
        tracer.samples = samples;
        tracer.stratified = stratified;
        tracer.color = color;
        tracer.phong = phong;
        tracer.reflections = reflections;
        tracer.shadows = shadows;
        tracer.maxraydepth = maxraydepth;
        tracer.verbose = verbose;
        tracer.createImage(argv[optind++]);
        std::cout << "Your image has been created!" << std::endl;
    }else{
        std::cout << "Error 404, try again later!" << std::endl;
    }



    return 0;
}
