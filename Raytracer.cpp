//William Whatley

#include "Raytracer.h"

Raytracer::Raytracer(Scene& scene){
    m_scene = scene;
}


//worldSpace makes it so that the camera is in a 3d plane
void Raytracer::worldSpace(){

    //all formulas from class
    vecW = m_scene.m_viewpoint->from - m_scene.m_viewpoint->at;
    vecW.normalize();

    vecU =  m_scene.m_viewpoint->up.cross(vecW);
    vecU.normalize();

    vecV = vecW.cross(vecU);
    vecV.normalize();

    Vector3d temp= m_scene.m_viewpoint->from - m_scene.m_viewpoint->at;
    magD = temp.norm();
}

void Raytracer::createImage(std::string output){
    //firstly call worldSpace to set up the vectors and magD
    worldSpace();

    //width and height of the image is based on the resolution
    //made constants so they aren't mistakenly changed
    const int WIDTH = m_scene.m_viewpoint->res[0];
    const int HEIGHT = m_scene.m_viewpoint->res[1];


    //field of view is the angle
    double fov = m_scene.m_viewpoint->angle;
    //converted into radians
    fov = fov * M_PI/double(180);
    //aspect ratio is width/height
    double aspect = double(WIDTH)/double(HEIGHT);

    //imagePlane and imageHeight set up
    const double planeWidth = 2 * magD * tan(fov / 2);
    const double planeHeight = planeWidth/ aspect;


    //pixel dimensions then calculated
    const double pixWidth = planeWidth / double(WIDTH);
    const double pixHeight = planeHeight / double(HEIGHT);

    //used to store the pixels
    unsigned char *pixels = new unsigned char [WIDTH * HEIGHT * 3];

    //image is then begun to be created
    for (int py = 0; py < HEIGHT; py++){
        for (int px = 0; px < WIDTH; px++){
            Vector3d color;
            //accounting for jittering
            if (stratified){
                //goes through to account for the samples
                for(int jy = 0; jy < samples; jy++){
                    for(int jx = 0; jx < samples; jx++){
                        //offsets
                        double ox = (jx + rand())/samples;
                        double oy = (jy + rand())/samples;

                        //ray direction calculated with jittering
                        double x = -planeWidth/2 + pixWidth * (px + ox);
                        double y = planeHeight/2 - pixHeight * (py + oy);

                        //direction set
                        Vector3d direction;
                        direction = (-1 * magD * vecW);
                        direction = direction + (x * vecU);
                        direction = direction + (y * vecV);
                        direction.normalize();

                        Ray ray(m_scene.m_viewpoint->from, direction);

                        color += colorSet(ray);
                    }
                }

                //color is then averaged and set
                color = color / (samples * samples);

                int index = (py * WIDTH + px) * 3;

                //set the colors
                double r = color[0];
                double g = color[1];
                double b = color[2];

                //now type cast and add to the array
                //in consecutive order to account for red blue and green values
                pixels[index] = static_cast<unsigned char>(std::min(255.0, std::max(0.0, r * 255)));
                pixels[index + 1] = static_cast<unsigned char>(std::min(255.0, std::max(0.0, g * 255)));
                pixels[index + 2] = static_cast<unsigned char>(std::min(255.0, std::max(0.0, b * 255)));
            }else{
                //direction declared
                Vector3d direction;

                //first need to calculate direction
                double x = -planeWidth/2 + pixWidth /2 + px * pixWidth;
                double y = planeHeight/2 - pixHeight /2 - py * pixHeight;

                //then create a vector based on that
                direction = (-1 * magD * vecW);
                direction = direction + (x * vecU);
                direction = direction + (y * vecV);
                direction = direction.normalized();

                Ray ray(m_scene.m_viewpoint->from, direction);

                //then determine if any object is hit using colorSet
                Vector3d color = colorSet(ray);

                //index is used to determine what pixels go where
                int index = (py * WIDTH + px) * 3;

                //set the colors
                double r = color[0];
                double g = color[1];
                double b = color[2];

                //now type cast and add to the array
                //in consecutive order to account for red blue and green values
                pixels[index] = static_cast<unsigned char>(std::min(255.0, std::max(0.0, r * 255)));
                pixels[index + 1] = static_cast<unsigned char>(std::min(255.0, std::max(0.0, g * 255)));
                pixels[index + 2] = static_cast<unsigned char>(std::min(255.0, std::max(0.0, b * 255)));
            }
        }
    }
    //makes it so that the new file is now printed
    const char* newFile = output.c_str();
    FILE *f = fopen(newFile,"wb");
    fprintf(f, "P6\n%d %d\n%d\n", WIDTH, HEIGHT, 255);
    fwrite(pixels, 1, HEIGHT*WIDTH*3, f);
    fclose(f);

    //pixels then deleted
    delete []pixels;
    
}

//colorSet determines if the color remains the background color or not
Vector3d Raytracer::colorSet(Ray ray){
    Surface* currSurface = nullptr;
    Hit hr;

    //checks to see if a shape is hit
    bool hit = isHit(ray, currSurface, hr);

    //if no hit occurs, return as it isn't useful
    if(!hit) return m_scene.m_background;

    hr.transfer(currSurface->m_fill);

    //local color is determined by the localLight
    Vector3d totalColor = localLight(ray, currSurface, hr);


    //if the m_fill[4] of the current surface (which corresponds to ks)
    //or raydepth is less than 5, then continue recursion
    if(hr.fill[4] > 0 && !(ray.depth > 5)){
        //eye and direction are set up
        Vector3d eye = (hr.inter + hr.norm * BIAS);
        Vector3d direction = ray.dir - 2 * (ray.dir.dot(hr.norm)) * hr.norm;
        direction.normalize();
        
        //reflectionRay then set up and used for recursion
        Ray reflectionRay(eye, direction, (ray.depth + 1));
        totalColor += hr.fill[4] * colorSet(reflectionRay);
    }

    return totalColor;

}

//localLight
//goes through the local light sources and determines the color vector
Vector3d Raytracer::localLight(Ray ray, Surface* currSurface, Hit& hit){
    //for loop goes through for each light source
    Vector3d localColor;
    Vector3d colorfill(hit.fill[0], hit.fill[1], hit.fill[2]);
    for(unsigned int i = 0; i < m_scene.m_lights.size(); i++){
        //the three vectors for Lambert shading are set up
        Vector3d l = (m_scene.m_lights[i]->coords - hit.inter).normalized();
        Vector3d v = -1 * ray.dir.normalized();
        Vector3d h = (l + v).normalized();

        //checks if the surface is a patch and that m_phong is being used
        if(phong && dynamic_cast<Patch*>(currSurface)){
            Patch* patch = static_cast<Patch*>(currSurface);
            Vector3d normal = patch->interpolate(hit);
            normal.normalize();
            hit.norm = normal;
        }   

        //shadow ray is created, and then is tested to see if the pt is in the shadow
        Ray shadowRay = Ray((hit.inter + hit.norm * BIAS), l);
        double distance = (m_scene.m_lights[i]->coords - hit.inter).norm();
        //if shadow test is false, then update
        if (!(shadowTest(shadowRay, distance))){
            //diffuse and specular are set
            double diffuse = (hit.norm.dot(l) > 0.0) ? hit.norm.dot(l) : 0.0;
            double max = (hit.norm.dot(h) > 0.0) ? hit.norm.dot(h) : 0.0;
            double specular = pow(max, hit.fill[5]);
            //intensity and colorfill set and then used to help calculate color
            double intensity = (1.0 / std::sqrt(m_scene.m_lights.size()));

            //update the colors
            localColor[0] += ((hit.fill[3] * colorfill[0] * diffuse) + (hit.fill[4] * specular)) * intensity;
            localColor[1] += ((hit.fill[3] * colorfill[1] * diffuse) + (hit.fill[4] * specular)) * intensity;
            localColor[2] += ((hit.fill[3] * colorfill[2] * diffuse) + (hit.fill[4] * specular)) * intensity;
        }
    }

    return localColor;
}

//shadowTest
//used to test if the shadow ray is represnted in the shadow
bool Raytracer::shadowTest(Ray ray, double distance){

    Hit hit; //hit set up but is not needed in this case

    //for loop goes through and determines if an intersection occurs
    for(unsigned int i = 0; i < m_scene.m_surfaces.size(); i++){
        //if an intersection occurs, then return true
        double bias = BIAS;
        if(m_scene.m_surfaces[i]->intersect(ray, bias, distance, hit)){
            return true;
        }
    }
    return false;
}

//isHit determines if an object is hit or not, and determines the type if so
bool Raytracer::isHit(Ray ray, Surface*& currSurface, Hit &hr){
    //originally delcares the shape hit as nothing
    bool state = false;

    //min set as 0 and max as infinity to start
    double min = 1e-6;
    double max = std::numeric_limits<double>::infinity();
    //goes through and sees if the ray hits any polygon
    for (const auto& surface : m_scene.m_surfaces){
        //if the verticies equals three, means we hit a triangle
        if(surface->intersect(ray, min, max, hr)){
            state = true;
            currSurface = surface;
            hr.raydepth = ray.depth;
            hr.view = ray.eye - hr.inter;
            hr.view.normalize();
        }
    }
    //if it goes through without hitting anything, returns NOTHING
    return state;
}
