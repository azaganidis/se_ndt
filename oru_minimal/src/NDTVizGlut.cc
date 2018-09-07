#include "ndt_visualisation/NDTVizGlut.hh"
// #include <assert.h>
// #include <unistd.h>
// #include <math.h>
// #include <iostream>

// // Needed for the wrapper calls.
static NDTVizGlut* glut3d_ptr = 0x0;

// Need some wrapper functions to handle the callbacks functions.
void win_reshape_(int w, int h) { glut3d_ptr->win_reshape (w,h); } 
void win_redraw_() { glut3d_ptr->win_redraw(); }
void win_key_(unsigned char key, int x, int y) { glut3d_ptr->win_key(key, x, y); }
void win_mouse_(int button, int state, int x, int y) { glut3d_ptr->win_mouse(button, state, x, y); }
void win_motion_(int x, int y) { glut3d_ptr->win_motion(x, y); }
void win_idle_() { glut3d_ptr->win_idle(); }  
void win_close_() { glut3d_ptr->win_close(); }

void * glthread(void * pParam)
{
    int argc=0;
    char** argv = NULL;
	      
    glutInit(&argc, argv);
	  
     // Create a window
     glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
     glutInitWindowSize(640,480);
     
     win = glutCreateWindow("NDTVizGlut");
     
     glEnable(GL_DEPTH_TEST);
     
     glEnable(GL_LIGHTING);
     glEnable(GL_LIGHT0);
     
     // Create light components
     GLfloat ambientLight[] = { 0.2f, 0.2f, 0.2f, 1.0f };
     GLfloat diffuseLight[] = { 0.8f, 0.8f, 0.8, 1.0f };
     GLfloat specularLight[] = { 0.5f, 0.5f, 0.5f, 1.0f };
     GLfloat position[] = { -1.5f, 1.0f, -4.0f, 1.0f };
     
     // Assign created components to GL_LIGHT0
     glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
     glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
     glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
     glLightfv(GL_LIGHT0, GL_POSITION, position);


     // enable color tracking
     glEnable(GL_COLOR_MATERIAL);
     // set material properties which will be assigned by glColor
     glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
     glClearColor(0.6f, 0.6f, 0.6f, 1.0f);

    glutReshapeFunc(win_reshape_);
     glutDisplayFunc(win_redraw_);
     glutKeyboardFunc(win_key_);
     glutMouseFunc(win_mouse_);
     glutMotionFunc(win_motion_);
     glutPassiveMotionFunc(NULL);

     // Idle loop callback
     glutIdleFunc(win_idle_);

     // Window close function
     glutCloseFunc(win_close_);
    /* Thread will loop here */

     glClearColor(0.0, 0.0, 0.0, 0.0);
     /* Thread will loop here */
     int f = 0;
     int n = glut3d_ptr->getMotionBlurFrames();
     if (glut3d_ptr->getFullScreen()) {
       glutFullScreen();
     }

#if 0
     while (true) {
         usleep(1000);
         for (int i = 0; i < 10; i++)
             glutMainLoopEvent();
         glut3d_ptr->update_cam();
         win_redraw_();
     }
     //glutMainLoop();
#endif

     while (true) {
       for (int i = 0; i < 2; i++) {
         glutMainLoopEvent();
       }
       
       glut3d_ptr->update_cam();
        while(glut3d_ptr->keyHit())
            usleep(10000);
       win_redraw_();
       if (n == 0) {
         glutSwapBuffers();
         usleep(10000);
       }
       else {
         // Performing motion bluring
         if(f == 0)
           glAccum(GL_LOAD, 1.0 / n);
         else
           glAccum(GL_ACCUM, 1.0 / n);
         f++;
         
         if(f >= n) {
           f = 0;
           glAccum(GL_RETURN, 1.0);
           glutSwapBuffers();
           usleep(10000);
         }
       }
     }

     return NULL;
}


NDTVizGlut::NDTVizGlut()
{
  aspect_ratio_factor = 1.;
  camera = &orbit_camera;

     // GUI settings
     gui_pause = 0;


     // cam_radius = 10.0f;
     // cam_azim = 0.5f;
     // cam_sweep_ang = 0.0f;

     // cam_sweep_origin.x = 0.0f;
     // cam_sweep_origin.y = 0.0f;
     // cam_sweep_origin.z = 0.0f;

     // cam_sweep_speed = 0.002;

     // cam_sweep = 0;

     glut3d_ptr = this;

     save_inc_counter = 0;
     do_save_inc = false;

     update_cam();

     open = true;
     full_screen = false;
     motion_blur_frames = 0;
}

NDTVizGlut::~NDTVizGlut()
{

}

void NDTVizGlut::setFullScreen(bool fs) 
{
     this->full_screen = fs;
}

bool NDTVizGlut::getFullScreen() const
{
     return this->full_screen;
}

void NDTVizGlut::setMotionBlurFrames(int f) 
{
  this->motion_blur_frames = f;
}

int NDTVizGlut::getMotionBlurFrames() const
{
  return this->motion_blur_frames;
}

void
NDTVizGlut::update_cam()
{
     glLoadIdentity();
     Eigen::Vector3f cp = camera->getPosition();
     Eigen::Vector3f fp = camera->getFocalPoint();
     Eigen::Vector3f up = camera->getUpVector();
     gluLookAt(cp[0], cp[1], cp[2],
	       fp[0], fp[1], fp[2],
	       up[0], up[1], up[2]);
}

void
NDTVizGlut::win_key(unsigned char key, int x, int y)
{
    this->pressed_keys.push_back(key);
}

// Mouse callback
void
NDTVizGlut::win_mouse(int button, int state, int x, int y)
{
//    std::cerr << "win_mouse - b:" << button << " s: " << state << "[" << x << "," << y << "]" << std::endl;
    camera->update_mouse(button, state, x, y);
    update_cam();
    //win_redraw();
    return;
}

void
NDTVizGlut::win_motion(int x, int y)
{
//    std::cerr << "win_motion : " << x << "," << y << std::endl;
    camera->update_motion(x, y);
    update_cam();
    //win_redraw();
    return;
}


// Handle window reshape events
void
NDTVizGlut::win_reshape(int width, int height)
{
     // Prevent a divide by zero, when window is too short
     // (you cant make a window of zero width).
     if(height == 0)
	  height = 1;
     
     float ratio = aspect_ratio_factor * width / height;
     // Reset the coordinate system before modifying
     glMatrixMode(GL_PROJECTION);
     glLoadIdentity();
     
     // Set the viewport to be the entire window
     glViewport(0, 0, width, height);
     
     // Set the clipping volume
     gluPerspective(45,ratio,1,1000);
     glMatrixMode(GL_MODELVIEW);
     
     update_cam();
     win_redraw();
     return;
}


// Redraw the window
void 
NDTVizGlut::win_redraw()
{
     glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

     this->draw_origin();

     // Do the drawing...
     this->draw();

     glutSwapBuffers();
     return;
}


void
NDTVizGlut::draw()
{
//     NDTVizGlut::draw();
     glColor3f(0.5, 0.5, 0.5);
     glBegin(GL_LINE_LOOP);
     glVertex3f(-2, -2, 0);
     glVertex3f(+2, -2, 0);
     glVertex3f(+2, +2, 0);
     glVertex3f(-2, +2, 0);
     glEnd();
     // Draw the objects.
     for (size_t i = 0; i < objects.size(); i++)
      {
          objects[i]->draw();
      }
}

// Idle callback
void 
NDTVizGlut::win_idle()
{
     if (!gui_pause)
     {
	  glutPostRedisplay();      
     }
     else
	  usleep(100000);
     
     return;
}

void
NDTVizGlut::win_close()
{
    std::cerr << "Window closed. " << std::endl;
    open = false;
}

bool 
NDTVizGlut::isOpen() const {
    return open;
}

bool
NDTVizGlut::keyHit() const {
    return !pressed_keys.empty();
}

void
NDTVizGlut::keyHitReset(){
    pressed_keys.pop_front();
}

unsigned char
NDTVizGlut::getPushedKey() {
    unsigned char ret = pressed_keys.front();
    return ret;
}

void *
start_glut_loop(void *ptr)
{
     if (ptr == 0x0)
     {
     }
     glutMainLoop();
     return 0x0;
}

void
NDTVizGlut::process_events()
{
//     if (do_save_inc)
//	  save_inc();
//     glutMainLoopEvent();
}

// Run the GUI
int 
NDTVizGlut::win_run(int *argc, char **argv)
{
    std::cerr << "win_run" << std::endl;
     // glutInit(argc, argv);
	  
     // // Create a window
     // glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
     // glutInitWindowSize(640,480);
     
     // win = glutCreateWindow("NDTVizGlut");
     
     // glEnable(GL_DEPTH_TEST);
     
     // glEnable(GL_LIGHTING);
     // glEnable(GL_LIGHT0);
     
     // // Create light components
     // GLfloat ambientLight[] = { 0.2f, 0.2f, 0.2f, 1.0f };
     // GLfloat diffuseLight[] = { 0.8f, 0.8f, 0.8, 1.0f };
     // GLfloat specularLight[] = { 0.5f, 0.5f, 0.5f, 1.0f };
     // GLfloat position[] = { -1.5f, 1.0f, -4.0f, 1.0f };
     
     // // Assign created components to GL_LIGHT0
     // glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
     // glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
     // glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
     // glLightfv(GL_LIGHT0, GL_POSITION, position);


     // // enable color tracking
     // glEnable(GL_COLOR_MATERIAL);
     // // set material properties which will be assigned by glColor
     // glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
     
//     glutIgnoreKeyRepeat(1);
     // glutReshapeFunc(win_reshape_);
     // glutDisplayFunc(win_redraw_);
     // glutKeyboardFunc(win_key_);
     // glutMouseFunc(win_mouse_);
     // glutMotionFunc(win_motion_);
     // glutPassiveMotionFunc(NULL);

     // // Idle loop callback
     // glutIdleFunc(win_idle_);

     // // Window close function
     // glutCloseFunc(win_close_);

     return 0;
}

void
NDTVizGlut::draw_origin()
{
     glBegin(GL_LINES);
     glColor3f(1, 0, 0);
     glVertex3f(0.0f, 0.0f, 0.0f); // origin of the FIRST line
     glVertex3f(1.0f, 0.0f, 0.0f); // ending point of the FIRST line
     glColor3f(0, 1, 0);
     glVertex3f(0.0f, 0.0f, 0.0f); // origin of the SECOND line
     glVertex3f(0.0f, 1.0f, 0.0f); // ending point of the SECOND line
     glColor3f(0, 0, 1);
     glVertex3f(0.0f, 0.0f, 0.0f);
     glVertex3f(0.0f, 0.0f, 1.0f);
     glEnd( ); 
}

int 
NDTVizGlut::save_inc()
{
    std::string file_name = "mov";// + getIntString(save_inc_counter, 4) + ".jpg";
     save_inc_counter++;
     return save(file_name);
}

int
NDTVizGlut::save(const std::string &fileName)
{
     // IplImage*          img;
     // GLint		viewport[4];
     // GLint		width, height;
     
     // glReadBuffer(GL_FRONT);
     // glGetIntegerv(GL_VIEWPORT, viewport);
     
     // width = viewport[2];
     // height = viewport[3];

     // // allocate the image
     // CvSize size = cvSize(width, height);
     // img = cvCreateImage( size, IPL_DEPTH_8U, 3 );
     
     // glFinish();
     // glPixelStorei(GL_PACK_ALIGNMENT, 4);
     // glPixelStorei(GL_PACK_ROW_LENGTH, 0);
     // glPixelStorei(GL_PACK_SKIP_ROWS, 0);
     // glPixelStorei(GL_PACK_SKIP_PIXELS, 0);
     
     
     // glReadPixels(0, 0, width, height,
     //    	  GL_BGR,
     //    	  GL_UNSIGNED_BYTE,
     //    	  img->imageData);

     // cvFlip(img); // The images are drawn upside down.
     // int ret_val = ocvSaveImage(img, fileName.c_str());
     // cvReleaseImage(&img);

     // return ret_val;
    return 1;
}

void
NDTVizGlut::repaint() {
    //this->win_redraw();

}

void
NDTVizGlut::clearScene() {
    objects.clear();
}

void
NDTVizGlut::setCameraPointingToPoint(double x, double y, double z) {
    camera->setFocalPoint(Eigen::Vector3f(x,y,z)); 
    update_cam();
}

void
NDTVizGlut::setCameraPosition(double x, double y, double z) {
  camera->setPosition(Eigen::Vector3f(x,y,z));
  update_cam();
}

void
NDTVizGlut::switchCamera(const std::string &type) {

  if (type == std::string("orbit")) {
    camera = &orbit_camera;
  }
  else if (type == std::string("fixed")) {
    camera = &fixed_camera;
  }
  else {
    assert(false);
  }
}

void
NDTVizGlut::setAspectRatioFactor(float factor) {
  aspect_ratio_factor = factor;
}

const NDTVizGlutCamera* 
NDTVizGlut::getCameraConstPtr() const {
  return camera;
}

NDTVizGlutCamera* 
NDTVizGlut::getCameraPtr() {
  return camera;
}
