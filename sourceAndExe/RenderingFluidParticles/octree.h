

struct WaterParticle;

struct OctreeNode {
  
  OctreeNode *parent;
  OctreeNode *child[8];
  
  double bbmin[3], bbmax[3];
  double inv_size[3];

  vector<double> surface_distance;

  OctreeNode(){
      parent = NULL;
      for(int i = 0; i < 8; i++) child[i] = NULL;
  }
};


class Octree {

    OctreeNode root;
    
    void ClipRaySegment(double *ray_origin, double *ray_direction, double *ray_direction_inv, double *bbmin, double *bbmax,
                            double &t_near, double &t_far);
    bool CheckNode(OctreeNode *node, double *org, double *dir, double *inv_dir, double t_near, double t_far, 
                   double *intersection, double *gradient, bool isComputeIntersection);
    void ReleaseNodeMemory(OctreeNode *node);
    void GetChildBBRange(double *bbmin, double *bbmax, double *new_bbmin, double *new_bbmax, int id);
    void ConstructOctree(OctreeNode *node, vector<WaterParticle> &water_particles, int depth);
    void AddChild(OctreeNode *node, int id);

    void DrawWireCube(OctreeNode *node, int level);

public:
    void InitOctree(vector<WaterParticle> &water_particles);
    void DeleteOctree();
    bool FindIntersection(double *org, double *dir, double *inv_dir, double *intersection, double *gradient, bool isComputeIntersection);
    void RenderOctreeNode();
};


