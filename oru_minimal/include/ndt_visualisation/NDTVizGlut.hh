#ifndef NDTVIZGLUT_HH
#define NDTVIZGLUT_HH

#include <GL/freeglut.h>
#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <deque>
#include <Eigen/Core>
#include <Eigen/Eigenvalues> 
#include <boost/thread.hpp>
 #include <assert.h>
 #include <unistd.h>
 #include <math.h>
 #include <iostream>

inline void checkOpenGLError()
{
  int	 openglErr;
  if ( ( openglErr= glGetError()) != GL_NO_ERROR )
    {
      const std::string sErr = std::string("OpenGL error: ") + std::string( (char*)gluErrorString(openglErr) );
      std::cerr << "[checkOpenGLError] " << sErr << std::endl;
    }
}

// typedef struct
// {
//      float x;
//      float y;
// } glut3d_vector2_t;

// typedef struct
// {
//      float x;
//      float y;
//      float z;
// } glut3d_vector3_t;


static int win;

//static void * glthread( void * pParam );
void * glthread( void * pParam );


static int create_gl_thread( void ) {
    
  pthread_t thread1;
     
  int iRet = pthread_create( &thread1, NULL, glthread, NULL );
     
  return iRet;
}


class NDTVizGlutColorPoint {
public:
  NDTVizGlutColorPoint() { }
  NDTVizGlutColorPoint(float x, float y, float z, float R, float G, float B) {
    pos[0] = x; pos[1] = y; pos[2] = z;
    col[0] = R; col[1] = G; col[2] = B;
  }
  GLfloat pos[3];
  GLfloat col[3];
};

class NDTVizGlutLine {
public:
  NDTVizGlutLine() { }
  NDTVizGlutLine(float x1, float y1, float z1, float x2, float y2, float z2) {
    pos1[0] = x1; pos1[1] = y1; pos1[2] = z1;
    pos2[0] = x2; pos2[1] = y2; pos2[2] = z2;
  }
  GLfloat pos1[3];
  GLfloat pos2[3];
};

class NDTVizGlutColor4f {
public:
  GLfloat R,G,B,A;
};

//! Interface class for all drawable objects
class NDTVizGlutObject
{
public:
  //! Draw
  virtual void draw() = 0;
};

                
class  NDTVizGlutText: public NDTVizGlutObject
{
public:
    float occupancy;

  NDTVizGlutText(float oc):occupancy(oc){};
  void draw(){
        glColor3f( 1, 1, 1 );
        glRasterPos2f(50, 50);
        int len, i;
        char OCCSTR[64];
        snprintf(OCCSTR,64, "occupancy: %f\0",occupancy);
        len = (int)strlen(OCCSTR);
        for (i = 0; i < len; i++) 
          glutBitmapCharacter(GLUT_BITMAP_9_BY_15, OCCSTR[i]);
  }
};
class NDTVizGlutPointCloudColor : public NDTVizGlutObject
{
public:
  NDTVizGlutPointCloudColor() {
    //m_pointSize = -1;
  }
  void draw() {
    //glEnable(GL_POINT_SMOOTH);
    //if (m_pointSize > 0.)
    //  glPointSize( m_pointSize );
    glBegin(GL_POINTS);
    for (size_t i = 0; i < pts.size(); i++) {
      glColor3fv(pts[i].col);
      glVertex3fv(pts[i].pos);
    }
    glEnd();
	checkOpenGLError();
  }
  void setPointSize(float size) { m_pointSize = size; }
  void clear() { pts.clear(); }
  void push_back(float x, float y, float z, float R, float G, float B) {
    pts.push_back(NDTVizGlutColorPoint(x,y,z,R,G,B));
  }
  NDTVizGlutColorPoint& getPoint(size_t i) { return pts[i]; }
private:
  float m_pointSize;
  std::vector<NDTVizGlutColorPoint> pts;
};

class NDTVizGlutSetOfLines : public NDTVizGlutObject
{
public:
  NDTVizGlutSetOfLines() {
    setColor4(1.0, 0.4, 0.0, 0.8);
    m_antiAliasing = true;
  }
  void draw() {
	glPushAttrib( GL_COLOR_BUFFER_BIT | GL_LINE_BIT );
	if (m_antiAliasing)
      glEnable(GL_LINE_SMOOTH);
    //	glLineWidth(m_lineWidth);
	checkOpenGLError();

	glDisable(GL_LIGHTING);  // Disable lights when drawing lines
	glBegin(GL_LINES);
	glColor4f(m_color.R,m_color.G,m_color.B,m_color.A);
	for (size_t i = 0; i < lines.size(); i++)
      {
        glVertex3fv(lines[i].pos1);
		glVertex3fv(lines[i].pos2);
      }
	glEnd();
	checkOpenGLError();
	glEnable(GL_LIGHTING);  // Disable lights when drawing lines

	// End of antialiasing:
    glPopAttrib();
  }
  void setPointSize(float size) { m_pointSize = size; }
  void clear() { lines.clear(); }
  void push_back(float x1, float y1, float z1, float x2, float y2, float z2) {
    lines.push_back(NDTVizGlutLine(x1, y1, z1, x2, y2, z2));
  }
  void appendLine(float x1, float y1, float z1, float x2, float y2, float z2) {
    this->push_back(x1, y1, z1, x2, y2, z2);
  }
  NDTVizGlutLine& getLine(size_t i) { return lines[i]; }
  void setColor(float R, float G, float B) { 
    m_color.R = R;
    m_color.G = G;
    m_color.B = B;
  }
  void setColor4(float R, float G, float B, float A) {
    m_color.R = R; m_color.G = G; m_color.B = B; m_color.A = A;
  }
  void setThickness(float thickness) {
    m_thickness = thickness;
  }

private:
  float m_lineWidth;
  NDTVizGlutColor4f m_color;
  bool m_antiAliasing;
  float m_pointSize;
  float m_thickness;
  std::vector<NDTVizGlutLine> lines;
};



class NDTVizGlutEllipsoid : public NDTVizGlutObject
{
public:
  NDTVizGlutEllipsoid() {
    m_quantiles = 3;
    m_lineWidth = 1;;
    m_drawSolid3D = true;
    m_3D_segments = 20;
    this->setColor4(0.5, 1.0, 0.1, 0.8);
  }
  NDTVizGlutEllipsoid(float quantiles, float lineWidth, bool drawSolid, int segments) {
    m_quantiles = quantiles;
    m_lineWidth = lineWidth;
    m_drawSolid3D = drawSolid;
    m_3D_segments = segments;
    this->setColor4(0.5, 1.0, 0.1, 0.8);
  }
  ~NDTVizGlutEllipsoid() { }
  void setLocation(double x, double y, double z) { 
    Eigen::Vector3d p(x,y,z);
    this->setPos(p);
  }
  void enableDrawSolid3D(bool s) { m_drawSolid3D = s; }
  void setPos(const Eigen::Vector3d &pos) { m_mean = pos; }
  void setCovMatrix(const Eigen::Matrix3d &cov) { this->setCov(cov); }
  void setCov(const Eigen::Matrix3d &cov) {
    m_cov = cov;
    const double d=m_cov.determinant();
	if (d==0 || d!=d) // Note: "d!=d" is a great test for invalid numbers, don't remove!
      {
        // All zeros:
        m_prevComputedCov = m_cov;
        m_eigVec = Eigen::Matrix3d::Zero(3,3);
        m_eigVal = Eigen::Matrix3d::Zero(3,3);
      }
	else
      {
        // Not null matrix: compute the eigen-vectors & values:
        m_prevComputedCov = m_cov;
        Eigen::EigenSolver<Eigen::Matrix3d> es(m_cov);
        //m_eigVal = es.eigenvalues();
        //m_eigVec = es.eigenvectors();

        m_eigVal = es.pseudoEigenvalueMatrix();
        m_eigVec = es.pseudoEigenvectors();
            
        m_eigVal = m_eigVal.cwiseSqrt();
        // Do the scale at render to avoid recomputing the m_eigVal for different m_quantiles
      }

  }

  void setColor4(float R, float G, float B, float A) {
    m_color.R = R; m_color.G = G; m_color.B = B; m_color.A = A;
  }

  void setColor(float R, float G, float B, float A) { setColor4(R, G, B, A); }

  void draw() {
    glPushMatrix();
    glTranslatef(m_mean(0), m_mean(1), m_mean(2));

    const size_t dim = m_cov.cols();
    if(m_eigVal(0,0) != 0.0 && m_eigVal(1,1) != 0.0 && (dim==2 || m_eigVal(2,2) != 0.0) && m_quantiles!=0.0)
      {
		glEnable(GL_BLEND);
		checkOpenGLError();
        glColor4f(m_color.R, m_color.G, m_color.B, m_color.A);
                
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		checkOpenGLError();
        glLineWidth(m_lineWidth);
		checkOpenGLError();

		if (dim==2)
          {
			// // ---------------------
			// //     2D ellipse
			// // ---------------------

			// /* Equivalent MATLAB code: 
			//  *
			//  * q=1;
			//  * [vec val]=eig(C); 
			//  * M=(q*val*vec)';  
			//  * R=M*[x;y];
			//  * xx=R(1,:);yy=R(2,:);  
			//  * plot(xx,yy), axis equal;
			//  */

			// double			ang;
			// unsigned int	i;

			// // Compute the new vectors for the ellipsoid:
            // Eigen::Matrix3d	M;
			// M.noalias() = double(m_quantiles) * m_eigVal * m_eigVec.adjoint();

			// glBegin( GL_LINE_LOOP );

			// // Compute the points of the 2D ellipse:
			// for (i=0,ang=0;i<m_2D_segments;i++,ang+= (M_2PI/m_2D_segments))
			// {
			// 	double ccos = cos(ang);
			// 	double ssin = sin(ang);

			// 	const float x = ccos * M.get_unsafe(0,0) + ssin * M.get_unsafe(1,0);
			// 	const float y = ccos * M.get_unsafe(0,1) + ssin * M.get_unsafe(1,1);

			// 	glVertex2f( x,y );
			// } // end for points on ellipse

			// glEnd();
          }
		else
          {
			// ---------------------
			//    3D ellipsoid
			// ---------------------
			GLfloat		mat[16];

			//  A homogeneous transformation matrix, in this order:
			//
			//     0  4  8  12
			//     1  5  9  13
			//     2  6  10 14
			//     3  7  11 15
			//
			mat[3] = mat[7] = mat[11] = 0;
			mat[15] = 1;
			mat[12] = mat[13] = mat[14] = 0;

			mat[0] = m_eigVec(0,0); mat[1] = m_eigVec(1,0); mat[2] = m_eigVec(2,0);	// New X-axis
			mat[4] = m_eigVec(0,1); mat[5] = m_eigVec(1,1); mat[6] = m_eigVec(2,1);	// New X-axis
			mat[8] = m_eigVec(0,2); mat[9] = m_eigVec(1,2); mat[10] = m_eigVec(2,2);	// New X-axis

			glDisable(GL_LIGHTING);
			//glEnable(GL_LIGHTING);
			//glEnable(GL_LIGHT0);
			glEnable(GL_COLOR_MATERIAL);
			glShadeModel(GL_SMOOTH);

			GLUquadricObj	*obj = gluNewQuadric();
			checkOpenGLError();

			gluQuadricDrawStyle( obj, m_drawSolid3D ? GLU_FILL : GLU_LINE);

			glPushMatrix();
			glMultMatrixf( mat );
            glScalef(m_eigVal(0,0)*m_quantiles,m_eigVal(1,1)*m_quantiles,m_eigVal(2,2)*m_quantiles);

			gluSphere( obj, 1,m_3D_segments,m_3D_segments);
			checkOpenGLError();

			glPopMatrix();

			gluDeleteQuadric(obj);
			checkOpenGLError();

			glDisable(GL_LIGHTING);
			glDisable(GL_LIGHT0);


          }


        glLineWidth(1.0f);
		glDisable(GL_BLEND);

      }
    glPopMatrix();
  }
private:
  Eigen::Matrix3d m_cov, m_prevComputedCov, m_eigVal, m_eigVec;
  Eigen::Vector3d m_mean;
  NDTVizGlutColor4f m_color;
  float m_quantiles;
  float m_lineWidth;
  bool m_drawSolid3D;
  int m_3D_segments;
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW  
};


class NDTVizGlutEllipsoids : public NDTVizGlutObject {
public:
  void draw() {
    for (size_t i = 0; i < objs.size(); i++) {
      objs[i].draw();
    }
  }
  void push_back(const NDTVizGlutEllipsoid &obj) { objs.push_back(obj); }
  void clear() { objs.clear(); objs.resize(0); }
private:
  std::vector<NDTVizGlutEllipsoid, Eigen::aligned_allocator<NDTVizGlutEllipsoid> > objs;
};

//! Interface class
class NDTVizGlutCamera {
public:
  virtual Eigen::Vector3f getPosition() const = 0;
  virtual Eigen::Vector3f getFocalPoint() const = 0;
  virtual Eigen::Vector3f getUpVector() const = 0;
  virtual void setFocalPoint(const Eigen::Vector3f &fp) = 0;
  virtual void setPosition(const Eigen::Vector3f &p) = 0;

  virtual void update_mouse(int button, int state, int x, int y) = 0;
  virtual void update_motion(int x, int y) = 0;
};

//! Implements a camera without (mouse)interaction
class NDTVizGlutFixedCamera : public NDTVizGlutCamera {
private:
  Eigen::Vector3f position_;
  Eigen::Vector3f focal_point_;
public:
  NDTVizGlutFixedCamera() {
    position_ = Eigen::Vector3f(0., 0., 1.);
    focal_point_ = Eigen::Vector3f(1., 0., 0.);
  }
  
  Eigen::Vector3f getPosition() const {
    return position_;
  }

  Eigen::Vector3f getFocalPoint() const {
    return focal_point_;
  }

  Eigen::Vector3f getUpVector() const {
    return Eigen::Vector3f(0., 0., 1.);
  }
  
  void setFocalPoint(const Eigen::Vector3f &fp) {
    focal_point_ = fp;
  }

  void setPosition(const Eigen::Vector3f &p) {
    position_ = p;
  }
  
  void update_mouse(int button, int state, int x, int y) {
    // Nothing
  }
  void update_motion(int x, int y) {
    // Nothing
  }
};

//! This implements a XY-orbit camera movement
class NDTVizGlutXYOrbitCamera : public NDTVizGlutCamera {
private: 
  int last_button_;
  int last_state_;
  Eigen::Vector2i last_point_;
  Eigen::Vector3f focal_point_;
  float distance;
  float yaw;
  float pitch;
public:
  NDTVizGlutXYOrbitCamera() {
    last_button_ = -1; 
    last_state_ = -1;
    distance = 30.;
    yaw = 1.;
    pitch = 0.7;
    focal_point_ = Eigen::Vector3f(0., 0., 0.);
  }

  Eigen::Vector3f getPosition() const {
    Eigen::Vector3f p = Eigen::Vector3f(distance*sin(pitch)*cos(yaw),
                                        distance*sin(pitch)*sin(yaw),
                                        distance*cos(pitch));
    return p + getFocalPoint();
  }
  Eigen::Vector3f getFocalPoint() const {
    return focal_point_;
  }
  void setFocalPoint(const Eigen::Vector3f &fp) {
    focal_point_ = fp;
  }
  void setPosition(const Eigen::Vector3f &p) {
    Eigen::Vector3f diff = p-getFocalPoint();
    distance = diff.norm();
    yaw = atan2(diff(1),diff(0));
    pitch = acos(diff(2) / distance);
    
    while (pitch > M_PI/2.) {
      pitch -= M_PI;
    }
    while (pitch < -M_PI/2.) {
      pitch += M_PI;
    }
      
    // std::cout << "distance : " << distance << std::endl;
    // std::cout << "yaw : " << yaw << std::endl;
    // std::cout << "pitch : " << pitch << std::endl;
  }
  Eigen::Vector3f getUpVector() const {
    return Eigen::Vector3f(0., 0., 1.);
  }
  float getPitchAngle() const { return pitch; }
  void setPitchAngle(double p) { pitch = p; } 
  void update_mouse(int button, int state, int x, int y) {
    last_button_ = button;
    last_state_ = state;
    last_point_[0] = x;
    last_point_[1] = y;

    // Wheel reports as button 3(scroll up) and button 4(scroll down)
    if ((button == 3) || (button == 4)) // It's a wheel event
      {
        // Each wheel event reports like a button click, GLUT_DOWN then GLUT_UP
        if (state != GLUT_UP) { // Disregard redundant GLUT_UP events
          if (button == 3) 
            distance *= 0.9;
          else
            distance *= 1.1;
        }
      }       
        
  }
  void update_motion(int x, int y) {
    // Left button in x -> change the yaw
    // Middle button -> move the position in x-dir,y-dir.
    // Rigth button in y -> move the position along z-dir.
    float dx = (x - last_point_[0])*0.01;
    float dy = (y - last_point_[1])*0.01;
    switch(last_button_) {
    case GLUT_LEFT_BUTTON:
      yaw += -dx;
      pitch += dy;
      if (pitch > M_PI/2.)
        pitch = M_PI/2.- 0.0001;
      if (pitch < -M_PI/2.)
        pitch = -M_PI/2.+ 0.0001;
      break;
    case GLUT_MIDDLE_BUTTON:
      dx *= distance*0.15;
      dy *= distance*0.15;
      focal_point_[1] += -cos(-yaw) * dx -sin(yaw) * dy;
      focal_point_[0] += -sin(-yaw) * dx -cos(yaw) * dy;
      //focal_point_[1] -= cos(yaw) * dx; // - sin(yaw) * dy;
      break;
    case GLUT_RIGHT_BUTTON:
      if (dy > 0)
        distance *= 1.05;
      else
        distance *= 0.95;
      break;
    default:
      break;
    };
    last_point_[0] = x;
    last_point_[1] = y;
  }
};


//! An OpenGL class to draw 3D stuff.
/*!
 * Based on the glut library. NOTE, this requires freeglut, not
 * the ordinary GLUT, since it requires glutMainLoopEvent() function.
 */
class NDTVizGlut {
public:
  //! Constructor
  NDTVizGlut();
  //! Destructor
  virtual ~NDTVizGlut();
  //! Run the GUI
  virtual int win_run(int *argc, char **argv);
  //! Key callback function
  virtual void win_key(unsigned char key, int x, int y);
  //! Mouse callback
  virtual void win_mouse(int button, int state, int x, int y);
  //! Mouse callback
  virtual void win_motion(int x, int y);
  //! Reshape events
  virtual void win_reshape(int width, int height);
  //! Redraw the window
  virtual void win_redraw();
  //! Idle callback
  virtual void win_idle();
  //! Close windo callback
  virtual void win_close();
  //! Process events
  virtual void process_events();

  virtual void start_main_loop() {
    glutMainLoop();
  }

  virtual void start_main_loop_own_thread() {
    // boost::thread workerThread(workerFunc);
    create_gl_thread();
    // if(pthread_create(&glut_event_processing_thread, NULL, ndt_viz_event_loop_thread, this)) {
	    
    //     std::cerr << "Error creating thread\n";
    // }
  }


  virtual void draw_origin();

  void setFullScreen(bool fs);
  bool getFullScreen() const;
  void setMotionBlurFrames(int f);
  int getMotionBlurFrames() const;

  //! Save an image (screenshot) of current view.
  int save(const std::string &fileName);

  //! Start/stop saving incrementally (to create movies).
  void set_save_inc_flag(bool flag) { do_save_inc = flag; }

  //! Add an object to draw.
  void addObject(NDTVizGlutObject* object) {
    // Only add objects if they are NOT in the objects.
    if ( std::find(objects.begin(), objects.end(), object)==objects.end() ) {
      objects.push_back(object);
    }
  }
     
  //! Set the drawing style
  //     void setDrawingStyle();

  //! 
  void repaint();

  void clearScene();

  void setCameraPointingToPoint(double x, double y, double z);
  void setCameraPosition(double x, double y, double z); //Not used?

  bool isOpen() const;

  bool keyHit() const;
  void keyHitReset();
  unsigned char getPushedKey();

  void update_cam();

  void switchCamera(const std::string &type);
  void setAspectRatioFactor(float ratio);

  const NDTVizGlutCamera* getCameraConstPtr() const;
  NDTVizGlutCamera* getCameraPtr();

protected:
  //! Put the code to draw here.
  /*!
   * This is called from the win_redraw function, which's also
   * draws the origin (0,0) of the 2D space.
   */
  virtual void draw();

  void cam_rotate();

  //! Saves a set of images.
  int save_inc();


  pthread_t glut_event_processing_thread;

  // GUI settings
  //     int win;
  int gui_pause;
  //     int viewport_width;
  //     int viewport_height;
  int show_samples;
  int show_grid;

  double start_time;

  Eigen::Vector3f cam_pos;
  Eigen::Vector3f cam_dir;

  // float cam_radius;
  // float cam_sweep_ang;
  // float cam_azim;
     
  // float cam_sweep_speed;

  // bool cam_sweep;

  int save_inc_counter;
  bool do_save_inc;

  bool open;

  // // Objects to draw
  std::vector<NDTVizGlutObject*> objects;
  NDTVizGlutCamera* camera;
  NDTVizGlutXYOrbitCamera orbit_camera;
  NDTVizGlutFixedCamera fixed_camera;

  std::deque<unsigned char> pressed_keys;

  float aspect_ratio_factor;
  bool full_screen;
  int motion_blur_frames;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////

// // // Needed for the wrapper calls.
// static NDTVizGlut* glut3d_ptr = 0x0;

// // Need some wrapper functions to handle the callbacks functions.
// void win_reshape_(int w, int h) { glut3d_ptr->win_reshape (w,h); } 
// void win_redraw_() { glut3d_ptr->win_redraw(); }
// void win_key_(unsigned char key, int x, int y) { glut3d_ptr->win_key(key, x, y); }
// void win_mouse_(int button, int state, int x, int y) { glut3d_ptr->win_mouse(button, state, x, y); }
// void win_motion_(int x, int y) { glut3d_ptr->win_motion(x, y); }
// void win_idle_() { glut3d_ptr->win_idle(); }  
// void win_close_() { glut3d_ptr->win_close(); }

// void * glthread(void * pParam)
// {
//   int argc=0;
//   char** argv = NULL;
	      
//   glutInit(&argc, argv);
	  
//   // Create a window
//   glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
//   glutInitWindowSize(640,480);
     
//   win = glutCreateWindow("NDTVizGlut");
     
//   glEnable(GL_DEPTH_TEST);
     
//   glEnable(GL_LIGHTING);
//   glEnable(GL_LIGHT0);
     
//   // Create light components
//   GLfloat ambientLight[] = { 0.2f, 0.2f, 0.2f, 1.0f };
//   GLfloat diffuseLight[] = { 0.8f, 0.8f, 0.8, 1.0f };
//   GLfloat specularLight[] = { 0.5f, 0.5f, 0.5f, 1.0f };
//   GLfloat position[] = { -1.5f, 1.0f, -4.0f, 1.0f };
     
//   // Assign created components to GL_LIGHT0
//   glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
//   glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
//   glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
//   glLightfv(GL_LIGHT0, GL_POSITION, position);


//   // enable color tracking
//   glEnable(GL_COLOR_MATERIAL);
//   // set material properties which will be assigned by glColor
//   glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);


//   glutReshapeFunc(win_reshape_);
//   glutDisplayFunc(win_redraw_);
//   glutKeyboardFunc(win_key_);
//   glutMouseFunc(win_mouse_);
//   glutMotionFunc(win_motion_);
//   glutPassiveMotionFunc(NULL);

//   // Idle loop callback
//   glutIdleFunc(win_idle_);

//   // Window close function
//   glutCloseFunc(win_close_);
//   /* Thread will loop here */

//   while (true) {
//     usleep(1000);
//     for (int i = 0; i < 10; i++)
//       glutMainLoopEvent();
//     glut3d_ptr->update_cam();
//     win_redraw_();
//   }
//   //glutMainLoop();
     
//   return NULL;
// }

#endif
