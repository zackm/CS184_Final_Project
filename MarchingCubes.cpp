//#include "MarchingCubes.h"
//
//
//MarchingCubes::MarchingCubes(void)
//{
//}
//
//
//MarchingCubes::~MarchingCubes(void)
//{
//}
//
///************************
//* Marching Cube Methods *
//************************/
///*
//Returns a list of pointers to the different triangles that need to be made for the marching cube case 
//corresponding to the input int. Should call some global list of all possible triangles (up to transformation)
//to make a particular triangle.
//*/
//vector<Triangle*> MarchingCubes::make_triangles(int i, int j, int k){
//	//cube define by the corner (i,j), (i,j+1), (i+1,j), and (i+1,j+1).
//
//	int cube_case = GRID_BOOL[i][j][k]*8+GRID_BOOL[i][j+1][k]*4+GRID_BOOL[i+1][j][k]*2+GRID_BOOL[i+1][j+1][k];
//
//	Vec3 vertex_1 = VERTEX_MATRIX[i][j][k];
//	Vec3 vertex_2 = VERTEX_MATRIX[i][j+1][k];
//	Vec3 vertex_3 = VERTEX_MATRIX[i+1][j][k];
//	Vec3 vertex_4 = VERTEX_MATRIX[i+1][j+1][k];
//	float weight_1 = GRID_DENSITY[i][j][k];
//	float weight_2 = GRID_DENSITY[i][j+1][k];
//	float weight_3 = GRID_DENSITY[i+1][j][k];
//	float weight_4 = GRID_DENSITY[i+1][j+1][k];
//
//	vector<Triangle*> tri_list;
//	Triangle* tri1, *tri2, *tri3;
//	Vec3 midpoint_1, midpoint_2;
//
//	if(cube_case==0){
//		//no triangles, no corners activated
//	}else if(cube_case==1){
//		//only corner 11 is turned on.
//		tri1 = new Triangle(vertex_4,vertex_2,vertex_3,weight_4,weight_2,weight_3,DENSITY_TOL);
//		tri_list.push_back(tri1);
//	}else if(cube_case==2){
//		//only corner 10 is turned on
//		tri1 = new Triangle(vertex_3,vertex_4,vertex_1,weight_3,weight_4,weight_1,DENSITY_TOL);
//		tri_list.push_back(tri1);
//	}else if(cube_case==3){
//		//10 and 11
//		midpoint_1 = vertex_3+(vertex_1-vertex_3)*((DENSITY_TOL-weight_3)/(weight_1-weight_3));
//		tri1 = new Triangle(vertex_3,vertex_4,midpoint_1);
//
//		midpoint_2 = vertex_4+(vertex_2-vertex_4)*((DENSITY_TOL-weight_4)/(weight_2-weight_4));
//		tri2 = new Triangle(vertex_4,midpoint_1,midpoint_2);
//
//		tri_list.push_back(tri1);
//		tri_list.push_back(tri2);
//	}else if(cube_case==4){
//		//only corner 01 is turned on
//		tri1 = new Triangle(vertex_2,vertex_1,vertex_4,weight_2,weight_1,weight_4,DENSITY_TOL);
//		tri_list.push_back(tri1);
//	}else if(cube_case==5){
//		//01 and 11 turned on
//		midpoint_1 = vertex_2+(vertex_1-vertex_2)*((DENSITY_TOL-weight_2)/(weight_1-weight_2));
//		tri1 = new Triangle(vertex_4,vertex_2,midpoint_1);
//
//		midpoint_2 = vertex_4+(vertex_3-vertex_4)*((DENSITY_TOL-weight_4)/(weight_3-weight_4));
//		tri2 = new Triangle(vertex_4,midpoint_1,midpoint_2);
//
//		tri_list.push_back(tri1);
//		tri_list.push_back(tri2);
//	}else if(cube_case==6){
//		//only 10 and 01 turned on
//		midpoint_1 = vertex_2+(vertex_1-vertex_2)*((DENSITY_TOL-weight_2)/(weight_1-weight_2));
//		midpoint_2 = vertex_2+(vertex_4-vertex_2)*((DENSITY_TOL-weight_2)/(weight_4-weight_2));
//		tri1 = new Triangle(vertex_2,midpoint_1,midpoint_2);
//
//		midpoint_1 = vertex_3+(vertex_4-vertex_3)*((DENSITY_TOL-weight_3)/(weight_4-weight_3));
//		midpoint_2 = vertex_3+(vertex_1-vertex_3)*((DENSITY_TOL-weight_3)/(weight_1-weight_3));
//		tri2 = new Triangle(vertex_3,midpoint_1,midpoint_2);
//
//		tri_list.push_back(tri1);
//		tri_list.push_back(tri2);
//	}else if(cube_case==7){
//		//only 10, 01, and 11 turned on
//		midpoint_1 = vertex_2+(vertex_1-vertex_2)*((DENSITY_TOL-weight_2)/(weight_1-weight_2));
//		tri1 = new Triangle(vertex_4,vertex_2,midpoint_1);
//
//		midpoint_2 = vertex_3+(vertex_1-vertex_3)*((DENSITY_TOL-weight_3)/(weight_1-weight_3));
//		tri2 = new Triangle(vertex_4,midpoint_1,midpoint_2);
//
//		tri3 = new Triangle(vertex_4,midpoint_2,vertex_3);
//
//		tri_list.push_back(tri1);
//		tri_list.push_back(tri2);
//		tri_list.push_back(tri3);
//	}else if(cube_case==8){
//		//only corner 00 is turned on
//		tri1 = new Triangle(vertex_1,vertex_3,vertex_2,weight_1,weight_3,weight_2,DENSITY_TOL);
//		tri_list.push_back(tri1);
//	}else if(cube_case==9){
//		//only 00 and 11 turned on
//		midpoint_1 = vertex_1+(vertex_3-vertex_1)*((DENSITY_TOL-weight_1)/(weight_3-weight_1));
//		midpoint_2 = vertex_1+(vertex_2-vertex_1)*((DENSITY_TOL-weight_1)/(weight_2-weight_1));
//		tri1 = new Triangle(vertex_1,midpoint_1,midpoint_2);
//
//		midpoint_1 = vertex_4+(vertex_2-vertex_4)*((DENSITY_TOL-weight_4)/(weight_2-weight_4));
//		midpoint_2 = vertex_4+(vertex_3-vertex_4)*((DENSITY_TOL-weight_4)/(weight_3-weight_4));
//		tri2 = new Triangle(vertex_4,midpoint_1,midpoint_2);
//
//		tri_list.push_back(tri1);
//		tri_list.push_back(tri2);
//	}else if(cube_case==10){
//		//only 00 and 10 turned on
//		midpoint_1 = vertex_3+(vertex_4-vertex_3)*((DENSITY_TOL-weight_3)/(weight_4-weight_3));
//		tri1 = new Triangle(vertex_1,vertex_3,midpoint_1);
//
//		midpoint_2 = vertex_1+(vertex_2-vertex_1)*((DENSITY_TOL-weight_1)/(weight_2-weight_1));
//		tri2 = new Triangle(vertex_1,midpoint_1,midpoint_2);
//
//		tri_list.push_back(tri1);
//		tri_list.push_back(tri2);
//	}else if(cube_case==11){
//		//00, 10, 11
//		midpoint_1 = vertex_4+(vertex_2-vertex_4)*((DENSITY_TOL-weight_4)/(weight_2-weight_4));
//		tri1 = new Triangle(vertex_3,vertex_4,midpoint_1);
//
//		midpoint_2 = vertex_1+(vertex_2-vertex_1)*((DENSITY_TOL-weight_1)/(weight_2-weight_1));
//		tri2 = new Triangle(vertex_3,midpoint_1,midpoint_2);
//
//		tri3 = new Triangle(vertex_3,midpoint_2,vertex_1);
//
//		tri_list.push_back(tri1);
//		tri_list.push_back(tri2);
//		tri_list.push_back(tri3);
//	}else if(cube_case==12){
//		//00 and 01 turned on
//		midpoint_1 = vertex_2+(vertex_4-vertex_2)*((DENSITY_TOL-weight_2)/(weight_4-weight_2));
//		tri1 = new Triangle(vertex_1,vertex_2,midpoint_1);
//
//		midpoint_2 = vertex_1+(vertex_3-vertex_1)*((DENSITY_TOL-weight_1)/(weight_3-weight_1));
//		tri2 = new Triangle(vertex_1,midpoint_1,midpoint_2);
//
//		tri_list.push_back(tri1);
//		tri_list.push_back(tri2);
//	}else if(cube_case==13){
//		//00, 01, 11
//		midpoint_1 = vertex_4+(vertex_3-vertex_4)*((DENSITY_TOL-weight_4)/(weight_3-weight_4));
//		tri1 = new Triangle(vertex_2,midpoint_1,vertex_4);
//
//		midpoint_2 = vertex_1+(vertex_3-vertex_1)*((DENSITY_TOL-weight_1)/(weight_3-weight_1));
//		tri2 = new Triangle(vertex_2,midpoint_2,midpoint_1);
//
//		tri3 = new Triangle(vertex_2,vertex_1,midpoint_2);
//
//		tri_list.push_back(tri1);
//		tri_list.push_back(tri2);
//		tri_list.push_back(tri3);
//	}else if(cube_case==14){
//		//00, 01, 10
//		midpoint_1 = vertex_2+(vertex_4-vertex_2)*((DENSITY_TOL-weight_2)/(weight_4-weight_2));
//		tri1 = new Triangle(vertex_1,midpoint_1,vertex_2);
//
//		midpoint_2 = vertex_3+(vertex_4-vertex_3)*((DENSITY_TOL-weight_3)/(weight_4-weight_3));
//		tri2 = new Triangle(vertex_1,midpoint_2,midpoint_1);
//
//		tri3 = new Triangle(vertex_1,vertex_3,midpoint_2);
//
//		tri_list.push_back(tri1);
//		tri_list.push_back(tri2);
//		tri_list.push_back(tri3);
//	}else if(cube_case==15){
//		tri1 = new Triangle(vertex_1,vertex_3,vertex_2);
//		tri2 = new Triangle(vertex_2,vertex_3,vertex_4);
//		tri_list.push_back(tri1);
//		tri_list.push_back(tri2);
//	}
//	return tri_list;
//}
//
///*
//Run marching cubes algorithm to generate triangles to render. Uses the particles in PARTICLES to
//calculate densities at corners of cubes (or squares in 2D). There are 16 cases for squares in
//the marching cubes algorithm, 256 for cubes (really 15 up to rotations and reflections).
//*/
//void MarchingCubes::marching_cubes(){
//	/*we need to break up screen into squares. We should store the squares in an efficient data
//	structure so we reuse values previously calculated through iterations.
//	*/
//
//	/*look here for cases http://users.polytech.unice.fr/~lingrand/MarchingCubes/algo.html
//	For now just making the 16 triangle types as is and enqueuing in list. Would be better to
//	make them dynamically.
//
//	Also, only doing uniform marching squares, not doing adaptive squares.
//	*/
//
//	VERTEX_MATRIX.clear();
//	GRID_DENSITY.clear();
//	GRID_BOOL.clear();
//	TRIANGLES.clear();
//
//	float error = .0001;
//	float step = CUBE_TOL;
//
//	//generate verticies and densities at those verticies
//	for (float x = CONTAINER.min.x; x<CONTAINER.max.x+error; x = x+step){
//		vector<vector<float> > yz_list;
//		vector<vector<Vec3> >vec_yz_list;
//		vector<vector<bool> > bool_yz_list;
//
//		for (float y = CONTAINER.min.y; y<CONTAINER.max.y+error; y = y+step){
//			vector<float> z_list;
//			vector<Vec3> vec_z_list;
//			vector<bool> bool_z_list;
//
//			for (float z = CONTAINER.min.z; z<CONTAINER.max.z+error; z = z+step){
//				Vec3 vertex(x,y,z);
//				float density = density_at_point(vertex);
//
//				z_list.push_back(density);
//				vec_z_list.push_back(vertex);
//				bool_z_list.push_back(density>DENSITY_TOL);
//			}
//			yz_list.push_back(z_list);
//			vec_yz_list.push_back(vec_z_list);
//			bool_yz_list.push_back(bool_z_list);
//		}
//
//		GRID_DENSITY.push_back(yz_list);
//		GRID_BOOL.push_back(bool_yz_list);
//		VERTEX_MATRIX.push_back(vec_yz_list);
//	}
//
//	//check marching cubes cases to add new triangles to list. It would be good to use bit operations for speed.
//	int m,n,p;//= GRID_DENSITY.size();//this implies uniform. Will need to change this.
//	m = n = p = GRID_BOOL.size();
//
//	int i,j,k; //i indicates the row, j indicates the column. 
//	i = j = k = 0;
//
//	bool squares_left = true;
//	vector<Triangle*> tri_list;
//	Vec3 vertex_1,vertex_2,vertex_3;
//	for (int i = 0; i<m-1;i++){
//		for (int j = 0; j<n-1;j++){
//			for (int k = 0; k<p-1;k++){
//				tri_list.clear();
//				tri_list = make_triangles(i,j,k);
//				TRIANGLES.insert(TRIANGLES.end(),tri_list.begin(),tri_list.end());
//
//				tri_list.clear();
//				tri_list = make_triangles(i,j,k+1);
//				TRIANGLES.insert(TRIANGLES.end(),tri_list.begin(),tri_list.end());
//			}	
//		}
//	}
//	//while (squares_left){ //need to change this later
//	//	tri_list.clear();
//	//	tri_list = make_triangles(i,j,k);
//	//	TRIANGLES.insert(TRIANGLES.end(),tri_list.begin(),tri_list.end());
//
//	//	tri_list.clear();
//	//	tri_list = make_triangles(i,j,k+1);
//	//	TRIANGLES.insert(TRIANGLES.end(),tri_list.begin(),tri_list.end());
//
//	//	//increment i and j and k or prepare to break loop
//	//	if(k==p-2){
//	//		if (j==n-2){
//	//			if (i==m-2){
//	//				squares_left = false;
//	//			}else{
//	//				j = 0;
//	//				i++;
//	//			}
//	//		}else{
//	//			j++;
//	//		}
//	//	}else{
//	//		k++;
//	//	}
//
//	//	if (j==n-2){
//	//		if (i==m-2){
//
//	//			squares_left = false;
//	//		}else{
//	//			j = 0;
//	//			i++;
//	//		}
//	//	}else{
//	//		j++;
//	//	}
//	//}
//}