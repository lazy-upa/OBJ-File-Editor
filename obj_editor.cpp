#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>

#include "raylib.h"

#define threshold 4
#define camConst 5

using namespace std;

struct coordinate {
	float p[4];
};

void scanLineFill(int x1, int y1, int x2, int y2, int x3, int y3, Color color) {
    int ymax = fmax(y1, fmax(y2, y3));
    int ymin = fmin(y1, fmin(y2, y3));

    for (int y = ymin; y <= ymax; y++) {
        int xLeft = 9999, xRight = -9999;

        if (y1 != y2 && y >= fmin(y1, y2) && y <= fmax(y1, y2))
            xLeft = xRight = x1 + (y - y1) * (x2 - x1) / (y2 - y1);
        if (y2 != y3 && y >= fmin(y2, y3) && y <= fmax(y2, y3)) {
            int x = x2 + (y - y2) * (x3 - x2) / (y3 - y2);
            xLeft = fmin(xLeft, x);
            xRight = fmax(xRight, x);
        }
        if (y3 != y1 && y >= fmin(y3, y1) && y <= fmax(y3, y1)) {
            int x = x3 + (y - y3) * (x1 - x3) / (y1 - y3);
            xLeft = fmin(xLeft, x);
            xRight = fmax(xRight, x);
        }

        for (int x = xLeft; x <= xRight; x++) {
            DrawPixel(x, y, color);
        }
    }
}

void line(int x1, int y1, int x2, int y2, Color color) {
    int dx = abs(x2 - x1), dy = abs(y2 - y1);
    int sx = (x1 < x2) ? 1 : -1;
    int sy = (y1 < y2) ? 1 : -1;
    int err = dx - dy;

    while (true) {
        DrawPixel(x1, y1, color);
        if (x1 == x2 && y1 == y2) break;
        int e2 = 2 * err;
        if (e2 > -dy) { err -= dy; x1 += sx; }
        if (e2 < dx) { err += dx; y1 += sy; }
    }
}

char temp = 'a';
int xmin, ymin, xmax, ymax, xc, yc;
bool editSel = false;
bool viewSel = false;
bool wireSel = false;

float l, r, t, b, n, f, aspect, fov_y;

class obj;
class matrix;
char projectionType = 'p';
char stateProjectionType ='p';

class normVector {
    public:
    float x;
    float y;
    float z;

    normVector(){
        x = 0.0;
        y = 0.0;
        z = 0.0;
    }

    normVector(float xx, float yy, float zz)
    {
        x = xx;
        y = yy;
        z = zz;
    }

    float magnitude(float x, float y, float z)
    {
        return sqrt(x*x + y*y + z*z);
    }

    void normalize()
    {

        float mag = magnitude(x, y, z);
        x = x/mag;
        y = y/mag;
        z = z/mag;
	}

    //cross product
    normVector operator * (normVector vec)
    {
        normVector res;

        res.x = y * vec.z - z * vec.y;
        res.y = z * vec.x - x * vec.z;
        res.z = x * vec.y - y * vec.x;

        return res;
    }
};

class matrix {
    private:
    float data[16];

    public:
    matrix(char c = 'i', float x = -1, float y = -1, float z = -1)
    {
        for(int i = 0; i < 16; i++)
        {
            data[i] = 0.0;
        }

        identity();

        switch(c)
        {
            case 't':
                translate(x,y,z);
                break;

            case 's':
                scale(x,y,z);
                break;

            case 'r':
                rotate((int)x, y);
                break;

            default:
                break;
        }
    }

    matrix(normVector u, normVector v, normVector n)
    {
    	identity();
        rotate(u,v,n);
    }

    matrix operator *(matrix other) const
    {
        matrix result;

        for (int row = 0; row < 4; row++){
            for (int col = 0; col < 4; col++) {
                    result.data[row * 4 + col] =
                    data[row * 4 + 0] * other.data[0 * 4 + col] +
                    data[row * 4 + 1] * other.data[1 * 4 + col] +
                    data[row * 4 + 2] * other.data[2 * 4 + col] +
                    data[row * 4 + 3] * other.data[3 * 4 + col];
            }
        }
        return result;
    }

    coordinate operator *(coordinate coor) const
    {

        coordinate temp = {};

        for(long int j = 0; j < 4; j++)
        {
            temp.p[j] = data[4 * j + 0] * coor.p[0] +
                data[4 * j + 1] * coor.p[1] +
                data[4 * j + 2] * coor.p[2] +
                data[4 * j + 3] * coor.p[3];
        }
		return temp;
    }

    matrix transpose() const
    {
        matrix trn;

        for(int i = 0; i < 4; i++)
        {
            trn.data[4 * i + 0] = data[0+ i];
            trn.data[4 * i + 1] = data[4+ i];
            trn.data[4 * i + 2] = data[8+ i];
            trn.data[4 * i + 3] = data[12+ i];
        }
        return trn;
    }

    void identity()
    {
        for(int i = 0; i < 16; i++)
        {
            data[i] = 0.0;
        }

        for(int i = 0; i < 4; i++)
        {
            data[4 * i + i] = 1;
        }
    }

    void translate(float tx = 0.0, float ty = 0.0, float tz = 0.0)
    {
        data[3] = tx;
        data[7] = ty;
        data[11] = tz;
    }

    void scale(float sx = 0.0, float sy = 0.0, float sz = 0.0)
    {
        data[0] = sx;
        data[5] = sy;
        data[10] = sz;
    }

    void rotate(int i, float theta)
    {
        theta *= 3.1415/180;

        float sine = sin(theta);
        float cosine = cos(theta);

        //-1 = rotation about x axis
        if(i == -1){
            data[5] = cosine;
            data[6] = -sine;

            data[9] = sine;
            data[10] = cosine;
        }
        //-2 = rotation about y axis
        if(i == -2){
            data[0] = cosine;
            data[2] = sine;

            data[8] = -sine;
            data[10] = cosine;
        }
        //-3 = rotation about z axis
        if(i == -3){
            data[0] = cosine;
            data[1] = -sine;

            data[4] = sine;
            data[5] = cosine;
        }

    }

    void rotate(normVector u, normVector v, normVector n)
    {
        data[0] = u.x;
        data[1] = u.y;
        data[2] = u.z;

        data[4] = v.x;
        data[5] = v.y;
        data[6] = v.z;

        data[8] = n.x;
        data[9] = n.y;
        data[10] = n.z;
    }

    friend class obj;
};

class obj{
    private:
    vector<float> vertices[4];
    vector<float> vn[4];
    vector <float> viewVertices[4];
    vector<float> pVertices[4];
    vector<float> pVn[4];
    vector<int> index;
    vector<int> faces[3];
    vector<int> faceNormal;
    vector<int> visible_face_normals;
    vector<float> selectedVertices[4];
    vector <int> visibleFaces[3];
    vector <float> midpointZ;

    public:
    float averageForCamera;

    obj() {
        averageForCamera = 0.0;
    }

    void transformNormals(matrix transform, char type = 'z')
    {
        float tNorm[4];

        for(unsigned long int i = 0, v = vn[0].size(); i < v; i++)
        {
            for(int j = 0; j < 4; j++)
            {
                tNorm[j] = transform.data[4 * j + 0] * vn[0][i] +
                    transform.data[4 * j + 1] * vn[1][i] +
                    transform.data[4 * j + 2] * vn[2][i] +
                    transform.data[4 * j + 3] * vn[3][i];
            }

            for(long int j = 0; j < 4; j++)
            {
                vn[j][i] = tNorm[j];
            }
        }
	}

    void transformSelectedVertex(matrix transform)
	{

        float temp[4];

        for(unsigned long int i = 0, v = selectedVertices[0].size(); i < v; i++)
        {
            for(long int j = 0; j < 4; j++)
            {
                temp[j] = transform.data[4 * j + 0] * selectedVertices[0][i] +
                    transform.data[4 * j + 1] * selectedVertices[1][i] +
                    transform.data[4 * j + 2] * selectedVertices[2][i] +
                    transform.data[4 * j + 3] * selectedVertices[3][i];
            }

            for(long int j = 0; j < 4; j++)
       		{
        		vertices[j][index[i]] = temp[j];
                selectedVertices[j][i] = temp[j];
        	}

    	}
	}

    void applyTransform(matrix transform, char type = 'z')
    {
        float temp[4];

        for(unsigned long int i = 0, v = vertices[0].size(); i < v; i++)
        {
            for(int j = 0; j < 4; j++)
            {
                temp[j] = transform.data[4 * j + 0] * vertices[0][i] +
                    transform.data[4 * j + 1] * vertices[1][i] +
                    transform.data[4 * j + 2] * vertices[2][i] +
                    transform.data[4 * j + 3] * vertices[3][i];
            }

            for(long int j = 0; j < 4; j++)
            {
                vertices[j][i] = temp[j];
            }
        }
    }

    void applyProjection(coordinate a, char type) {
        float at[3] = {0, 0, 0};

        normVector w(a.p[0] - at[0], a.p[1] - at[1], a.p[2] - at[2]);
        w.normalize();

        normVector u;
        if (abs(w.y) > 0.99) {
            u =  normVector(0,0,-1)* w;
            u.normalize();
        }
        else {
            u =  normVector(0,1,0)* w;
            u.normalize();
        }

        normVector v = w * u;
        v.normalize();

        matrix rot(u, v, w);
        matrix camera = rot * matrix('t', -a.p[0], -a.p[1], -a.p[2]);

        float tnorm[4], temp[4];
        for(unsigned long int i = 0, vSize = vn[0].size(); i < vSize; i++)
        {
            for(int j = 0; j < 4; j++)
            {
                tnorm[j] = rot.data[4 * j + 0] * vn[0][i] +
                    rot.data[4 * j + 1] * vn[1][i] +
                    rot.data[4 * j + 2] * vn[2][i] +
                    rot.data[4 * j + 3] * vn[3][i];

                pVn[j][i] = tnorm[j];
            }

        }

        for (unsigned long int i = 0, vSize = vertices[0].size(); i < vSize; i++) {
            for (int j = 0; j < 4; j++) {
                temp[j] = camera.data[4 * j + 0] * vertices[0][i] +
                        camera.data[4 * j + 1] * vertices[1][i] +
                        camera.data[4 * j + 2] * vertices[2][i] +
                        camera.data[4 * j + 3] * vertices[3][i];
                viewVertices[j][i] = temp[j];
            }

        }

        if (type == 'o'){
            l = -averageForCamera * camConst;
            r = -l;
            t = -l;
            b = -t;
            n = 1;
            f = 1000;

        matrix viewToClip;

        viewToClip.data[0] = 2 / (r - l);
        viewToClip.data[3] = -((r + l) / (r - l));
        viewToClip.data[5] = 2 / (t - b);
        viewToClip.data[7] = -((t + b) / (t - b));
        viewToClip.data[10] = -2 / (f - n); // Changed sign from -2 to 2
        viewToClip.data[11] = -(f + n) / (f - n); // Removed negative sign

            matrix ctoclip = viewToClip * camera;

            for (unsigned long int i = 0, vSize = vertices[0].size(); i < vSize; i++) {
                for (int j = 0; j < 4; j++) {
                    temp[j] = ctoclip.data[4 * j + 0] * vertices[0][i] +
                            ctoclip.data[4 * j + 1] * vertices[1][i] +
                            ctoclip.data[4 * j + 2] * vertices[2][i] +
                            ctoclip.data[4 * j + 3] * vertices[3][i];
                    pVertices[j][i] = temp[j];
                }
            }




            for (unsigned long int i = 0, vSize = vertices[0].size(); i < vSize; i++) {
                float x_ndc = pVertices[0][i];
                float y_ndc = pVertices[1][i];
                float z_ndc = pVertices[2][i];



                // if (x_ndc < -1) x_ndc = -1;
                // if (x_ndc > 1) x_ndc = 1;
                // if (y_ndc < -1) y_ndc = -1;
                // if (y_ndc > 1) y_ndc = 1;
                // if (z_ndc < -1) z_ndc = -1;
                // if (z_ndc > 1) z_ndc = 1;

                pVertices[0][i] = (float)(xmax - 0) / 2 * x_ndc + (float)(xmax + xmin) / 2;
                pVertices[1][i] = (float)(ymax - ymin) / 2 * (-y_ndc) + (float)(ymax + ymin) / 2;
                pVertices[2][i] = (f - n) / 2 * z_ndc + (f + n) / 2;


            }

        }
        if (type == 'p') {
            n = 1.0f;
            f = 1000.0f;

            fov_y = 60.0f * (3.14159265f / 180.0f);
            t =  n * tan(fov_y / 2.0f);
            b = -t;

            aspect = 1.0f;
            r = t * aspect;
            l = -r;
            matrix viewToClip;
            viewToClip.data[0] = (2 * n) / (r - l);
            viewToClip.data[2] = ((r + l) / (r - l));
            viewToClip.data[5] = 2 * n / (t - b);
            viewToClip.data[6] = ((t + b) / (t - b));
            viewToClip.data[11] = -2 * f *n / (f - n);
            viewToClip.data[10] = -((f + n) / (f - n));
            viewToClip.data[14] = -1;
            viewToClip.data[15] = 0;

            matrix ctoclip = viewToClip * camera;

            for (unsigned long int i = 0, vSize = vertices[0].size(); i < vSize; i++) {
                for (int j = 0; j < 4; j++) {
                    temp[j] = ctoclip.data[4 * j + 0] * vertices[0][i] +
                            ctoclip.data[4 * j + 1] * vertices[1][i] +
                            ctoclip.data[4 * j + 2] * vertices[2][i] +
                            ctoclip.data[4 * j + 3] * vertices[3][i];
                    pVertices[j][i] = temp[j];
                }

                float w_coord = pVertices[3][i];
                if (w_coord != 0) {
                pVertices[0][i] /= w_coord;
                pVertices[1][i] /= w_coord;
                pVertices[2][i] /= w_coord;
            }
            }




            for (unsigned long int i = 0, vSize = vertices[0].size(); i < vSize; i++) {
                float x_ndc = pVertices[0][i];
                float y_ndc = pVertices[1][i];
                float z_ndc = pVertices[2][i];



                if (x_ndc < -1) x_ndc = -1;
                if (x_ndc > 1) x_ndc = 1;
                if (y_ndc < -1) y_ndc = -1;
                if (y_ndc > 1) y_ndc = 1;
                if (z_ndc < -1) z_ndc = -1;
                if (z_ndc > 1) z_ndc = 1;

                pVertices[0][i] = (float)(xmax - 0) / 2 * x_ndc + (float)(xmax + xmin) / 2;
                pVertices[1][i] = (float)(ymax - ymin) / 2 * (-y_ndc) + (float)(ymax + ymin) / 2;
                pVertices[2][i] = (f - n) / 2 * z_ndc + (f + n) / 2;


            }

        }
    }

    void drawObject(coordinate c, coordinate at)
    {
        applyProjection(c, projectionType);
        // display

        visibleFaces[0].clear();
        visibleFaces[1].clear();
        visibleFaces[2].clear();
        midpointZ.clear();
        visible_face_normals.clear();

        for(unsigned long int i = 0, f = faces[0].size(); i < f; i++)
        {

            if(!wireSel){
                if (projectionType == 'p') {
                    normVector d(viewVertices[0][faces[0][i]], viewVertices[1][faces[0][i]] , viewVertices[2][faces[0][i]]);
                    d.normalize();

                    float dot = d.x * pVn[0][faceNormal[i]] + d.y * pVn[1][faceNormal[i]] + d.z * pVn[2][faceNormal[i]];

                    if(dot  > 0)
                    {
                        continue;
                    }
                }
                else if (projectionType == 'o') {
                    if (pVn[2][faceNormal[i]] <= 0) {
                        continue;
                    }
                }
            }

            visibleFaces[0].push_back(faces[0][i]);
            visibleFaces[1].push_back(faces[1][i]);
            visibleFaces[2].push_back(faces[2][i]);
            midpointZ.push_back((pVertices[2][faces[0][i]] + pVertices[2][faces[1][i]] +pVertices[2][faces[2][i]] )/3);
            visible_face_normals.push_back(faceNormal[i]);

        }

        sortVisibleFaces();
        // cout << visibleFaces[0].size() << endl;
        // cout << faces[0].size() << endl;

        for(unsigned long int i = 0, f = visibleFaces[0].size(); i < f; i++)
        {
            float xf1 = pVertices[0][visibleFaces[0][i]];
            float yf1 = pVertices[1][visibleFaces[0][i]];
            float xf2 = pVertices[0][visibleFaces[1][i]];
            float yf2 = pVertices[1][visibleFaces[1][i]];
            float xf3 = pVertices[0][visibleFaces[2][i]];
            float yf3 = pVertices[1][visibleFaces[2][i]];

            if (xf1 > xmax && xf2 > xmax && xf3 > xmax) {
                if (yf1 > ymax && yf2 > ymax && yf3 > xmax) {
                    continue;
                }
                if (yf1 < ymin && yf2 < ymin && yf3 < ymin) {
                    continue;
                }
            }

            if (xf1 < xmin && xf2 < xmin && xf3 < xmin) {
                if (yf1 > ymax && yf2 > ymax && yf3 > xmax) {
                    continue;
                }
                if (yf1 < ymin && yf2 < ymin && yf3 < ymin) {
                    continue;
                }
            }

            // at least one of the points are outside but some are inside
            if (xf1 > xmax) xf1=xmax;
            if (yf1 > ymax) yf1=ymax;
            if (xf1 < xmin) xf1=xmin;
            if (yf1 < ymin) yf1=ymin;

            if (xf2 > xmax) xf2=xmax;
            if (yf2 > ymax) yf2=ymax;
            if (xf2 < xmin) xf2=xmin;
            if (yf2 < ymin) yf2=ymin;

            if (xf2 > xmax) xf2=xmax;
            if (yf2 > ymax) yf2=ymax;
            if (xf2 < xmin) xf2=xmin;
            if (yf2 < ymin) yf2=ymin;

            if (xf3 > xmax) xf3=xmax;
            if (yf3 > ymax) yf3=ymax;
            if (xf3 < xmin) xf3=xmin;
            if (yf3 < ymin) yf3=ymin;


            if(!wireSel)
            {
                int cl = 40 + 180 * abs(pVn[2][visible_face_normals[i]]);

                 Color fillColor = {static_cast<unsigned char> (cl), static_cast<unsigned char> (cl), static_cast<unsigned char> (cl), 255};

                 scanLineFill(xf1, yf1, xf2, yf2, xf3, yf3, fillColor);
            }

            if(wireSel && !editSel)
            {
                line(xf1, yf1, xf2, yf2, WHITE);
                line(xf2, yf2, xf3, yf3, WHITE);
                line(xf3, yf3, xf1, yf1, WHITE);
            }

            if(!wireSel && editSel){
                line(xf1, yf1, xf2, yf2, BLACK);
                line(xf2, yf2, xf3, yf3, BLACK);
                line(xf3, yf3, xf1, yf1, BLACK);
                dotsOnVertex(visibleFaces[0][i]);
                dotsOnVertex(visibleFaces[1][i]);
                dotsOnVertex(visibleFaces[2][i]);
            }

            if (wireSel && editSel) {
                line(xf1, yf1, xf2, yf2, WHITE);
                line(xf2, yf2, xf3, yf3, WHITE);
                line(xf3, yf3, xf1, yf1, WHITE);
                dotsOnVertex(visibleFaces[0][i], WHITE);
                dotsOnVertex(visibleFaces[1][i], WHITE);
                dotsOnVertex(visibleFaces[2][i], WHITE);
            }

        }

        if (editSel) {
            for(unsigned long int i = 0, v = selectedVertices[0].size(); i < v; i++)
            {
                dotsOnVertex(i, Color({255,0,0,255}), 's');
            }
        }


    }

    void sortVisibleFaces()
    {
        float tempZ;
        int tempFace[3];
        int tempNorm;

        unsigned long int  v = visibleFaces[0].size();
        for(unsigned long int i = 0; i < v; i++)
        {
            for(unsigned long int j = i+1; j < v; j++)
            {
                if(abs(midpointZ[i]) < abs(midpointZ[j]))
                {
                    //swap
                    tempZ = midpointZ[i];
                    midpointZ[i] = midpointZ[j];
                    midpointZ[j] = tempZ;

                    tempNorm = visible_face_normals[i];
                    visible_face_normals[i] = visible_face_normals[j];
                    visible_face_normals[j] = tempNorm;

                    tempFace[0] = visibleFaces[0][i];
                    tempFace[1] = visibleFaces[1][i];
                    tempFace[2] = visibleFaces[2][i];

                    visibleFaces[0][i] = visibleFaces[0][j];
                    visibleFaces[1][i] = visibleFaces[1][j];
                    visibleFaces[2][i] = visibleFaces[2][j];

                    visibleFaces[0][j] = tempFace[0];
                    visibleFaces[1][j] = tempFace[1];
                    visibleFaces[2][j] = tempFace[2];
                }
            }
        }
    }

    void loadObject(const string &address){
        vertices[0].clear();
        vertices[1].clear();
        vertices[2].clear();
        vertices[3].clear();
        faces[0].clear();
        faces[1].clear();
        faces[2].clear();
        averageForCamera = 0;
		std::ifstream f(address.c_str());
		if(!f){
			std::cerr << "Error loading the obj file!";
		}

		std::string line;
		while(std::getline(f, line)){

			if(line.empty() || line[0] == '#'){
				continue;
			}

	        if(line[0] == 'v' && line[1] == ' '){
	            float x, y, z;
	            sscanf(line.c_str(), "v %f %f %f", &x, &y, &z);
                averageForCamera += (abs(x)+abs(y)+abs(z));
	            vertices[0].push_back(x);
	            vertices[1].push_back(y);
	            vertices[2].push_back(z);
	            vertices[3].push_back(1.0f);
	        }




	        if(line[0] == 'f' && line[1] == ' '){
	            int a, b, c;
	            int na, nb, nc;
	            int ta, tb, tc;

	            if(line.find("//") != std::string::npos){
                    temp = 'd';
	                sscanf(line.c_str(), "f %d//%d %d//%d %d//%d", &a, &na, &b, &nb, &c, &nc);
	                faces[0].push_back(a-1);
	                faces[1].push_back(b-1);
	                faces[2].push_back(c-1);

	            }

	            else if(line.find("/") != std::string::npos){
                    temp = 's';
	                sscanf(line.c_str(), "f %d/%d/%d %d/%d/%d %d/%d/%d", &a, &ta, &na, &b, &tb, &nb, &c, &tc, &nc);
	                faces[0].push_back(a-1);
	                faces[1].push_back(b-1);
	                faces[2].push_back(c-1);

	            }

	            else{
                    temp = 'n';
	                sscanf(line.c_str(), "f %d %d %d", &a, &b, &c);
	                faces[0].push_back(a-1);
	                faces[1].push_back(b-1);
	                faces[2].push_back(c-1);
	            }
	        }
		}
        averageForCamera /= (3 * vertices[0].size());
		for (long int j = 0; j < 4; j++)
		{
    		pVertices[j].resize(vertices[j].size());
            viewVertices[j].resize(vertices[j].size());
            pVn[j].resize(vn[j].size());
		}

		cout << "File Loaded Successfully!" << endl;
		f.close();
        calculateNormals();
	}

         void saveObject(const std::string &filename) {
            std::ofstream o(filename.c_str());

            if (!o) {
                std::cerr << "Error opening file for writing: " << filename << std::endl;
                return;
            }

            for (size_t i = 0; i < vertices[0].size(); i++) {
                o << "v " << vertices[0][i] << " " << vertices[1][i] << " " << vertices[2][i] << "\n";
            }

            o << std::endl;

            for (size_t i = 0; i < vn[0].size(); i++) {
                o << "vn " << vn[0][i] << " " << vn[1][i] << " " << vn[2][i] << "\n";
            }

            o << std::endl;

            for (size_t i = 0; i < faces[0].size(); i++) {
                    o << "f " << faces[0][i]+1 << "//" << faceNormal[i] +1<< " "
                      << faces[1][i]+1 << "//" << faceNormal[i] +1<< " "
                      << faces[2][i] +1<< "//" << faceNormal[i] +1<< "\n";
            }

            o.close();
            std::cout << "File saved successfully!" << std::endl;
        }

    void selectVertices(unsigned long int ind)
    {
        selectedVertices[0].push_back(vertices[0][ind]);
        selectedVertices[1].push_back(vertices[1][ind]);
        selectedVertices[2].push_back(vertices[2][ind]);
        selectedVertices[3].push_back(vertices[3][ind]);
        index.push_back(static_cast <int> (ind));
    }

    void clearSelectedVertices()
    {
        selectedVertices[0].clear();
        selectedVertices[1].clear();
        selectedVertices[2].clear();
        selectedVertices[3].clear();
        index.clear();
    }

    void dotsOnVertex(int i, Color clr = {0,0,0,255}, char c = 'p')
    {

        if(c == 's')
        {
            for(int j = -2; j < 3; j++)
            {
                for(int k = -2; k < 3; k++)
                {
                    DrawPixel((int)pVertices[0][index[i]]+j, (int)pVertices[1][index[i]]+k, clr);
                }
            }
            return;
        }

        for(int j = -2; j < 3; j++)
        {
            for(int k = -2; k < 3; k++)
            {
                DrawPixel((int)pVertices[0][i]+j, (int)pVertices[1][i]+k, clr);
            }
        }


    }

    bool closeToVertex(float xMouse, float yMouse, float xVertex, float yVertex) {
        float dx = xVertex - xMouse;
        float dy = yVertex - yMouse;
        return ((abs(dx) <= threshold) && (abs(dy) <= threshold));
    }

    void checking(float mx, float my) {
            for (unsigned long int i = 0, v = visibleFaces[0].size(); i < v; i++) {
                if (closeToVertex(mx, my, pVertices[0][visibleFaces[0][i]], pVertices[1][visibleFaces[0][i]])) {
                    selectVertices(visibleFaces[0][i]);
                    return;
                }
                if (closeToVertex(mx, my, pVertices[0][visibleFaces[1][i]], pVertices[1][visibleFaces[1][i]])) {
                    selectVertices(visibleFaces[1][i]);
                    return;
                }
                if (closeToVertex(mx, my, pVertices[0][visibleFaces[2][i]], pVertices[1][visibleFaces[2][i]])) {
                    selectVertices(visibleFaces[2][i]);
                    return;
                }
            }
    }

    void calculateNormals()
    {


        for(int k = 0; k < 4; k++)
        {
            vn[k].clear();
        }
        faceNormal.clear();

        for(unsigned long int i = 0, f = faces[0].size(); i < f; i++)
        {
            normVector a(vertices[0][faces[1][i]] - vertices[0][faces[0][i]], vertices[1][faces[1][i]] - vertices[1][faces[0][i]], vertices[2][faces[1][i]] - vertices[2][faces[0][i]]);

            normVector b(vertices[0][faces[2][i]] - vertices[0][faces[0][i]], vertices[1][faces[2][i]] - vertices[1][faces[0][i]], vertices[2][faces[2][i]] - vertices[2][faces[0][i]]);


            normVector n = a * b;
            n.normalize();

                vn[0].push_back(n.x);
                vn[1].push_back(n.y);
                vn[2].push_back(n.z);
                vn[3].push_back(1.0f);

            faceNormal.push_back((int)i);
        }


        for (int j = 0; j < 4; j++)
		{
            pVn[j].resize(vn[j].size());
		}
    }

    float averageX() {
        if (selectedVertices[0].size() == 0) return 0.0f;
        float sum = 0.0f;
        for (int i = 0; i < selectedVertices[0].size(); i++) {
            sum += selectedVertices[0][i];
        }
        return sum / selectedVertices[0].size();
    }

    float averageY() {
        if (selectedVertices[0].size() == 0) return 0.0f;
        float sum = 0.0f;
        for (int i = 0; i < selectedVertices[0].size(); i++) {
            sum += selectedVertices[1][i];
        }
        return sum / selectedVertices[0].size();
    }

    float averageZ() {
        if (selectedVertices[0].size() == 0) return 0.0f;
        float sum = 0.0f;
        for (int i = 0; i < selectedVertices[0].size(); i++) {
            sum += selectedVertices[2][i];
        }
        return sum / selectedVertices[0].size();
    }

    int getSelectedVerticesSize() {
        return selectedVertices[0].size();
    }
};

void DrawAxes(int vMode) {
    int viewWidth = xmax-xmin;
    int viewHeight = ymax-ymin;
    int centerX = (xmax + xmin) / 2;
    int centerY = (ymax + ymin) / 2;

    int axisLength = viewHeight / 2;



    if (vMode == 1) {
        line(centerX - axisLength, centerY, centerX + axisLength, centerY, RED);
        line(centerX, centerY - axisLength, centerX, centerY + axisLength, GREEN);

    }
    else if (vMode == 2) {
        line(centerX - axisLength, centerY, centerX + axisLength, centerY, RED);
        line(centerX, centerY - axisLength, centerX, centerY + axisLength, BLUE);

    }
    else if (vMode == 3) {
        line(centerX - axisLength, centerY, centerX + axisLength, centerY, GREEN);
        line(centerX, centerY - axisLength, centerX, centerY + axisLength, BLUE);

    }
}

int main(){
	// Initialize graphics window
	int sW = 1200;
	int sH = 960;


	xmax = 850;
	ymax = 930;

	xmin = 50;
	ymin = 130;

	xc = (xmax+xmin)/2;
	yc = (ymax+ymin)/2;

	coordinate c = {5,5,5,1};

    coordinate at = {};
	at.p[0] = 0;
	at.p[1] = 0;
	at.p[2] = 0;
	at.p[3] = 1;

    obj *ptr = nullptr;

	obj cube;
    cube.loadObject("cube.obj");
    //cube.applyTransform(matrix('s', 100, 100, 100));
    cube.applyTransform(matrix('t', -0.5, -0.5,-0.5));
    cube.calculateNormals();

    obj teapot;
    teapot.loadObject("teapot.obj");
    //teapot.applyTransform(matrix('s', 50, 50, 50));

    obj sphere;
    sphere.loadObject("sphere.obj");
    //sphere.applyTransform(matrix('s', 100, 100, 100));

    obj cylinder;
    cylinder.loadObject("cylinder.obj");
    obj custom;

    // c = {averageForCamera * 5, averageForCamera * 5 , averageForCamera * 5, 1};

    InitWindow(sW, sH, "3D-Window");

    SetTargetFPS(20);               // Set our game to run at 60 frames-per-second
    //--------------------------------------------------------------------------------------

        Rectangle buttonCube = {900, 100, 300, 100};
        Rectangle buttonTeapot = {900, 0, 300, 100};
        Rectangle buttonCylinder = {900, 300, 300, 100};
        Rectangle buttonSphere = {900, 200, 300, 100};
        Rectangle buttonCustom = {900, 400, 300, 100};
        Rectangle buttonExport = {900, 500, 300, 100};
        Rectangle textCustom = {900, 700, 300, 50};
        Rectangle textExport = {900, 700, 300, 50};
        //Rectangle clipWindow = {xmin, ymin, xmax-xmin, ymax-ymin};

        char customText[32] = "";
        int customTextCount = 0;
        bool customTextActive = false;

        char exportText[32] = "";
        int exportTextCount = 0;
        bool exportTextActive = false;

        Rectangle EditMode = {0, 0, 300, 100};
        Rectangle ViewMode = {300, 0, 300, 100};
        Rectangle WireMode = {600, 0, 300, 100};

        bool cubeSel = false;
        bool teapotSel = false;
        bool exportSel = false;
        bool cylinderSel = false;
        bool sphereSel = false;
        bool customSel = false;

        bool dragging = false;

        bool rotating = false;
        int rotationAxis = 0;

        bool translating = false;
        int translationAxis = 0;

        bool scaling = false;
        int scalingAxis = 0;

        float dxx = c.p[0] - at.p[0];
        float dyy = c.p[1] - at.p[1];
        float dzz = c.p[2] - at.p[2];
        float dd;

        int vMode = 0;

    // Main game loop
    while (!WindowShouldClose())    // Detect window close button or ESC key
    {

        Vector2 mousePoint = GetMousePosition();

        //buttonCubeHovered = CheckCollisionPointRec(mousePoint, buttonCube);
        //buttonOctHovered = CheckCollisionPointRec(mousePoint, buttonOct);

        if(CheckCollisionPointRec(mousePoint, buttonCustom) && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)){
            ClearBackground(BLACK);
            dragging = false;
            cubeSel = false;
            teapotSel = false;
            editSel = false;
            viewSel = false;
            wireSel = false;
            exportSel = false;
            cylinderSel = false;
            sphereSel = false;
            customSel = true;
            customTextActive = true;
            exportTextActive = false;
        }

        if(customTextActive){
            int key = GetCharPressed();
            while(key > 0){
                if((key >= 32 && key <= 125) && (customTextCount < 32)){
                    customText[customTextCount] = char(key);
                    customText[customTextCount+1] = '\0';
                    customTextCount++;
                }
                key = GetCharPressed();
            }

            if(IsKeyPressed(KEY_BACKSPACE)){
                customTextCount--;
                customText[customTextCount] = '\0';
            }

            if(IsKeyPressed(KEY_ENTER)){
                customTextActive = false;

                ptr = &custom;
                ptr -> loadObject(customText);
                c = {ptr -> averageForCamera * camConst, ptr -> averageForCamera * camConst, ptr -> averageForCamera * camConst, 1};

            }
        }

         if (CheckCollisionPointRec(mousePoint, buttonExport) && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
            ClearBackground(BLACK);
            dragging = false;
            cubeSel = false;
            teapotSel = false;
            editSel = false;
            viewSel = false;
            wireSel = false;
            exportSel = true;
            cylinderSel = false;
            sphereSel = false;
            customSel = false;
            exportTextActive = true; // Enable text input for export
            customTextActive = false;
        }

        if (exportTextActive) {
            int key = GetCharPressed();
            while (key > 0) {
                if ((key >= 32 && key <= 125) && (exportTextCount < 32)) {
                    exportText[exportTextCount] = char(key);
                    exportText[exportTextCount + 1] = '\0';
                    exportTextCount++;
                }
                key = GetCharPressed();
            }

            if (IsKeyPressed(KEY_BACKSPACE)) {
                if (exportTextCount > 0) {
                    exportTextCount--;
                    exportText[exportTextCount] = '\0';
                }
            }

            if (IsKeyPressed(KEY_ENTER)) {
                exportTextActive = false;
                if (ptr != nullptr){
                    ptr -> saveObject(exportText);
                    cout << "saved\n";
                }
                else{
                    //DrawText("Saving Failed!", 100, 100,20, BLACK);
                    std::cout << "Saving Failed!\n";
                //ptr -> applyTransform(matrix('s', 100, 100, 100));
                //ptr -> applyTransform(matrix('t', -50, -50,-50));
            }
        }
        }

        if(CheckCollisionPointRec(mousePoint, buttonCube) && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)){
            ptr = &cube;
            ClearBackground(BLACK);
            cubeSel = true;
            teapotSel= false;
            // editSel = false;
            // viewSel = false;
            // wireSel = false;
            exportSel = false;
            cylinderSel = false;
            sphereSel = false;
            customSel = false;
            exportTextActive = false;
            customTextActive = false;
            c  = {cube.averageForCamera * camConst, cube.averageForCamera * camConst, cube.averageForCamera * camConst, 1};

        }

        if(CheckCollisionPointRec(mousePoint, buttonTeapot) && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)){
            ptr = &teapot;
            ClearBackground(BLACK);
            cubeSel = false;
            teapotSel = true;
            // editSel = false;
            // viewSel = false;
            // wireSel = false;
            exportSel = false;
            cylinderSel = false;
            sphereSel = false;
            customSel = false;
            exportTextActive = false;
            customTextActive = false;
            c  = {teapot.averageForCamera * camConst, teapot.averageForCamera * camConst, teapot.averageForCamera * camConst, 1};
        }

        if(CheckCollisionPointRec(mousePoint, buttonCylinder) && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)){
            ptr = &cylinder;
            ClearBackground(BLACK);
            dragging = false;
            cubeSel = false;
            teapotSel = false;
            // editSel = false;
            // viewSel = false;
            // wireSel = false;
            exportSel = false;
            cylinderSel = true;
            sphereSel = false;
            customSel = false;
            exportTextActive = false;
            customTextActive = false;
            c  = {cylinder.averageForCamera * camConst, cylinder.averageForCamera * camConst, cylinder.averageForCamera * camConst, 1};
        }

        if(CheckCollisionPointRec(mousePoint, buttonSphere) && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)){
            ptr = &sphere;
            ClearBackground(BLACK);
            dragging = false;
            cubeSel = false;
            teapotSel = false;
            // editSel = false;
            // viewSel = false;
            // wireSel = false;
            exportSel = false;
            cylinderSel = false;
            sphereSel = true;
            customSel = false;
            exportTextActive = false;
            customTextActive = false;
            c  = {sphere.averageForCamera * camConst, sphere.averageForCamera * camConst, sphere.averageForCamera * camConst, 1};
        }

        if(CheckCollisionPointRec(mousePoint, EditMode) && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)){
            ClearBackground(BLACK);
            editSel = true;
            viewSel = false;
            // cubeSel = false;
            // teapotSel = false;
            // wireSel = false;
            exportSel = false;
            // cylinderSel = false;
            // sphereSel = false;
            // customSel = false;
            exportTextActive = false;
            customTextActive = false;
        }

        if(CheckCollisionPointRec(mousePoint, ViewMode) && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)){
            ClearBackground(BLACK);
            dragging = false;
            // cubeSel = false;
            // teapotSel = false;
            editSel = false;
            viewSel = true;
            wireSel = false;
            exportSel = false;
            // cylinderSel = false;
            // sphereSel = false;
            // customSel = false;
            exportTextActive = false;
            customTextActive = false;
        }

        if(CheckCollisionPointRec(mousePoint, WireMode) && IsMouseButtonPressed(MOUSE_LEFT_BUTTON)){
            ClearBackground(BLACK);
            dragging = false;
            // cubeSel = false;
            // teapotSel = false;
            // editSel = false;
            viewSel = false;
            wireSel = true;
            exportSel = false;
            // cylinderSel = false;
            // sphereSel = false;
            // customSel = false;
            exportTextActive = false;
            customTextActive = false;
        }

        if (IsMouseButtonDown(MOUSE_LEFT_BUTTON)) {
            if (!dragging) {
                dragging = true; // Start dragging
            }

            if (dragging) {
                Vector2 deltaMouse = GetMouseDelta();
                float dx = deltaMouse.x;
                float dy = deltaMouse.y;

                normVector w(c.p[0]-at.p[0], c.p[1]-at.p[1], c.p[2]-at.p[2]);
                w.normalize();
                normVector u, v;


                u = w * normVector(0,1,0);
                u.normalize();

                v = w * u;


                matrix rot(u,v,w);

                matrix cam = matrix('r', -1, dy) * rot;
                matrix inv_rot = rot.transpose();
                matrix final = inv_rot * cam;

                c = final * c;
                c = matrix('r', -2, -dx) * c;



                //ptr->applyProjection(cam); // Apply transformation to the selected object
            }
        } else {
            dragging = false; // Stop dragging when the button is not down
        }

        if(IsMouseButtonPressed(MOUSE_RIGHT_BUTTON) && editSel == true)
        {
            Vector2 mouse = GetMousePosition();
            if(ptr != nullptr) ptr -> checking((float)mouse.x, (float)mouse.y);
        }

        if(IsKeyPressed(KEY_Q))
        {
            ptr->clearSelectedVertices();
        }

        if (IsKeyPressed(KEY_ONE))
        {    vMode = 1;
            dd = sqrt(c.p[0]*c.p[0] + c.p[1]*c.p[1] + c.p[2]*c.p[2]);
            c.p[0] = at.p[0] + ptr->averageForCamera * camConst;
            c.p[1] = at.p[1];
            c.p[2] = at.p[2];
        }
        if (IsKeyPressed(KEY_TWO)){
            vMode = 2;
            dd = sqrt(c.p[0]*c.p[0] + c.p[1]*c.p[1] + c.p[2]*c.p[2]);
            c.p[0] = at.p[0];
            c.p[1] = at.p[1] + ptr->averageForCamera * camConst;
            c.p[2] = at.p[2];
        }
        if (IsKeyPressed(KEY_THREE)){ vMode = 3;
            dd = sqrt(c.p[0]*c.p[0] + c.p[1]*c.p[1] + c.p[2]*c.p[2]);
            c.p[0] = at.p[0];
            c.p[1] = at.p[1];
            c.p[2] = at.p[2] + ptr->averageForCamera * camConst;
        }

        if (vMode != 0) {
            DrawAxes(vMode);
            projectionType ='o';
            if(IsMouseButtonPressed(MOUSE_LEFT_BUTTON))
            {
                if (vMode == 1) {c = {dd,0,0,1};}
                if (vMode == 2) {c = {0,dd,0,1};}
                if (vMode == 3) {c = {0,0,dd,1};}
                vMode=0;
                projectionType = stateProjectionType;

            }
        }

        if (IsKeyPressed(KEY_R))
        {
            dragging = false;
            rotating = true;
            rotationAxis = 0;
        }
        if (rotating) {
            if (IsKeyPressed(KEY_X))
                rotationAxis = 1;
            else if (IsKeyPressed(KEY_Y))
                rotationAxis = 2;
            else if (IsKeyPressed(KEY_Z))
                rotationAxis = 3;
        }
        if (rotating && rotationAxis != 0 && ptr != nullptr) {
            if(ptr->getSelectedVerticesSize() == 0 && editSel){
                Vector2 deltarot = GetMouseDelta();
                float rotationspeed = (deltarot.x + deltarot.y) / 2;
                matrix rotMatrix = matrix('r', -rotationAxis, rotationspeed);
                ptr->applyTransform(rotMatrix);
                ptr->transformNormals(rotMatrix);
            }
            else {
                Vector2 deltarotvert = GetMouseDelta();
                float rotationspeed = (deltarotvert.x + deltarotvert.y) / 2.0f;

                matrix translateToOrigin = matrix('t', -ptr->averageX(),-ptr->averageY(),-ptr->averageZ());
                matrix rotMatrix = matrix('r', -rotationAxis, rotationspeed);
                matrix translateBack = matrix('t', ptr->averageX(),ptr->averageY(),ptr->averageZ());

                matrix finalMatrix = translateBack * rotMatrix * translateToOrigin;

                ptr->transformSelectedVertex(finalMatrix);
                ptr->calculateNormals();
            }
        }

        if (IsKeyPressed(KEY_T))
        {
            dragging = false;
            translating = true;
            translationAxis = 0;
        }
        if (translating)
        {
            if (IsKeyPressed(KEY_X))
                translationAxis = 1;
            else if (IsKeyPressed(KEY_Y))
                translationAxis = 2;
            else if (IsKeyPressed(KEY_Z))
                translationAxis = 3;
        }
        if (translating && translationAxis != 0 && ptr != nullptr)
        {
            if(ptr->getSelectedVerticesSize() == 0 && editSel){
                Vector2 deltaTrans = GetMouseDelta();
                float transAmount = (translationAxis == 1) ? deltaTrans.x*0.1 :(translationAxis == 2) ? deltaTrans.y*0.1 :(deltaTrans.x*0.1 + deltaTrans.y*0.1) / 2;
                float tx = (translationAxis == 1) ? transAmount : 0.0f;
                float ty = (translationAxis == 2) ? transAmount : 0.0f;
                float tz = (translationAxis == 3) ? transAmount : 0.0f;
                matrix transMatrix = matrix('t', tx, ty, tz);
                ptr->applyTransform(transMatrix);
            }
            else{
                Vector2 deltaTransvert = GetMouseDelta();
                float transAmount = (translationAxis==1)? deltaTransvert.x*0.1: (translationAxis ==2)?
                deltaTransvert.y*0.1 :(deltaTransvert.x*0.1 + deltaTransvert.y*0.1)/2;
                float tx = (translationAxis ==1)?transAmount:0.0f;
                float  ty = (translationAxis == 2) ?transAmount :0.0f;
                float tz = (translationAxis ==3)? transAmount :0.0f;
                matrix transMatrixvert = matrix('t', tx, ty, tz);
                ptr->transformSelectedVertex(transMatrixvert);
                ptr->calculateNormals();
            }
        }

       if(projectionType =='p'){
        float scrollDelta = GetMouseWheelMove();
        if (scrollDelta != 0.0f) {
            float zoomSpeed = 0.2f;


            float dx = c.p[0] - at.p[0];
            float dy = c.p[1] - at.p[1];
            float dz = c.p[2] - at.p[2];


            float distance = sqrt(dx * dx + dy * dy + dz * dz);


            float newDistance = distance - zoomSpeed * scrollDelta;
            if(newDistance < 1.0f) newDistance = 1.0f;


            normVector viewDir(dx, dy, dz);
            viewDir.normalize();

            c.p[0] = at.p[0] + viewDir.x * newDistance;
            c.p[1] = at.p[1] + viewDir.y * newDistance;
            c.p[2] = at.p[2] + viewDir.z * newDistance;

            }
        }

        if (IsKeyPressed(KEY_S))
        {
            dragging = false;
            scaling = true;
        }
        if (scaling)
        {
            if (IsKeyPressed(KEY_X))
                scalingAxis = 1;
            else if (IsKeyPressed(KEY_Y))
                scalingAxis = 2;
            else if (IsKeyPressed(KEY_Z))
                scalingAxis = 3;
        }

        if (scaling && ptr != nullptr)
        {
            //  && scalingAxis != 0
            Vector2 deltascale = GetMouseDelta();
            float scaleFactor = 1.0f + (deltascale.x + deltascale.y) * 0.01f;
            float sx, sy, sz;

            if (scalingAxis == 0) {
                sx = scaleFactor;
                sy = scaleFactor;
                sz = scaleFactor;
            }
            else {
                sx = (scalingAxis == 1) ? scaleFactor : 1.0f;
                sy = (scalingAxis == 2) ? scaleFactor : 1.0f;
                sz = (scalingAxis == 3) ? scaleFactor : 1.0f;
            }


            if (ptr->getSelectedVerticesSize() == 0 && editSel) {

                matrix scaleMatrix = matrix('s', sx, sy, sz);
                ptr->applyTransform(scaleMatrix);
                ptr->calculateNormals();
            }
            else {

                float avgX = ptr->averageX();
                float avgY = ptr->averageY();
                float avgZ = ptr->averageZ();

                matrix translateToOrigin = matrix('t', -avgX, -avgY, -avgZ);
                matrix scaleMatrix = matrix('s', sx, sy, sz);
                matrix translateBack = matrix('t', avgX, avgY, avgZ);

                matrix finalMatrix = translateBack * scaleMatrix * translateToOrigin;

                ptr->transformSelectedVertex(finalMatrix);
                ptr->calculateNormals();

            }

        }

        if (IsMouseButtonReleased(MOUSE_LEFT_BUTTON))
        {
            dragging = false;
            rotating = false;
            translating = false;
            scaling = false;
            rotationAxis = translationAxis = scalingAxis = 0;
            // ptr -> clearSelectedVertices();
        }
        if (IsKeyPressed(KEY_FIVE)) {
            if (projectionType == 'o'){
                projectionType = 'p';
                stateProjectionType = 'p';

            }
            else
            {
                projectionType = 'o';
                stateProjectionType='o';
            }
        }


        BeginDrawing();
        ClearBackground(BLACK);
            if(customTextActive) {
                //ClearBackground(BLACK);
                DrawRectangleRec(textCustom, LIGHTGRAY);
                DrawText("Enter obj file's name:", 900, 680, 20, WHITE);
                DrawText(customText, textCustom.x + 5, textCustom.y + 15, 20, BLACK);
            }

            if(exportTextActive) {
                //ClearBackground(BLACK);
                DrawRectangleRec(textExport, LIGHTGRAY);
                DrawText("Enter export obj file's name:", 900, 680, 20, WHITE);
                DrawText(exportText, textExport.x + 5, textExport.y + 15, 20, BLACK);
            }


            DrawRectangleRec(buttonCube, (cubeSel ? GREEN : LIGHTGRAY));
            DrawRectangleRec(buttonTeapot, (teapotSel ? GREEN : LIGHTGRAY));
            DrawRectangleRec(EditMode, (editSel ? GREEN : LIGHTGRAY));
            DrawRectangleRec(ViewMode, (viewSel ? GREEN : LIGHTGRAY));
            DrawRectangleRec(WireMode, (wireSel ? GREEN : LIGHTGRAY));
            DrawRectangleRec(buttonCustom, (customSel ? GREEN : LIGHTGRAY));
            DrawRectangleRec(buttonCylinder, (cylinderSel? GREEN : LIGHTGRAY));
            DrawRectangleRec(buttonSphere, (sphereSel ? GREEN : LIGHTGRAY));
            DrawRectangleRec(buttonExport, (exportSel? GREEN : LIGHTGRAY));
            //DrawRectangleRec(clipWindow, BLACK);
            DrawText("Cube", buttonCube.x + 120, buttonCube.y + 40, 30, BLACK);
            DrawText("Teapot", buttonTeapot.x + 120, buttonTeapot.y + 40, 30, BLACK);
            DrawText("Edit Mode", EditMode.x + 120, EditMode.y + 40, 20, BLACK);
            DrawText("View Mode", ViewMode.x + 120, ViewMode.y + 40, 20, BLACK);
            DrawText("Wire Mode", WireMode.x + 120, WireMode.y + 40, 20, BLACK);
            DrawText("Cylinder", buttonCylinder.x + 120, buttonCylinder.y + 40, 30, BLACK);
            DrawText("Sphere", buttonSphere.x + 120, buttonSphere.y + 40, 30, BLACK);
            DrawText("Custom", buttonCustom.x + 120, buttonCustom.y + 40, 30, BLACK);
            DrawText("Export", buttonExport.x + 120, buttonExport.y + 40, 30, BLACK);

            //For Border
            line(xmin, ymin, xmin, ymax, WHITE);
            line(xmax, ymin, xmax, ymax, WHITE);
            line(xmin, ymin, xmax, ymin, WHITE);
            line(xmin, ymax, xmax, ymax, WHITE);

            if(projectionType == 'o'){
                DrawText("ORTHOGRAPHIC", 1000, 930, 20, WHITE);
            }
            else if(projectionType == 'p'){
                DrawText("PERSPECTIVE", 1000, 930, 20, WHITE);
            }

            if (vMode == 1){
                DrawText("FRONT VIEW",700, 150, 20, WHITE);
            }
            else if (vMode == 2){
                DrawText("TOP VIEW",700, 150, 20, WHITE);
            }

            else if(vMode == 3){
                DrawText("SIDE VIEW",700, 150, 20, WHITE);
            }

            ClearBackground(Color ({25,25,25,255}));
            // c = matrix('r', -2, 1) * c;


            ClearBackground(Color ({25,25,25,255}));
            // c = matrix('r', -2, 1) * c;


            if(ptr != nullptr)
            {
                ptr->drawObject(c, at);
            }
            //cube.applyTransform(matrix('r', -2, 1));


        EndDrawing();

        }
    CloseWindow();
	return 0;
}