// Real-Time Physics Tutorials
// Brandon Pelfrey
// SPH Fluid Simulation
//
// g++ tutorial_3.cpp -lglut -O3 -fopenmp -o tutorial_3

#include <GL/glut.h>
#include <omp.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
using namespace std;

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/tee.hpp>
#include <boost/chrono.hpp>
#include <boost/chrono/chrono_io.hpp>
using namespace boost;

#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>
using namespace glm;

// --------------------------------------------------------------------
// A structure for holding two neighboring particles and their weighted distances
struct Neighbor
{
    int j;
    float q, q2;
};

// The Particle structure holding all of the relevant information.
struct Particle
{
    vec2 pos;
    float r, g, b;

    vec2 pos_old;
    vec2 vel;
    vec2 force;
    float mass;
    float rho;
    float rho_near;
    float press;
    float press_near;
    float sigma;
    float beta;
    vector< Neighbor > neighbors;
};

// --------------------------------------------------------------------
float G = .02f * .25;           // Gravitational Constant for our simulation

float spacing = 2.f;            // Spacing of particles
float k = spacing / 1000.0f;    // Far pressure weight
float k_near = k*10;            // Near pressure weight
float rest_density = 3;         // Rest Density
float r= spacing * 1.25f;       // Radius of Support
float rsq = r * r;              // ... squared for performance stuff

float SIM_W = 50;               // The size of the world
float bottom = 0;               // The floor of the world

// Our collection of particles
vector< Particle > particles;

// Mouse attractor
vec2 attractor(999,999);
bool attracting = false;


// Between [0,1]
float rand01()
{
    return (float)rand() * (1.f / RAND_MAX);
}

// Between [a,b]
float randab(float a, float b)
{
    return a + (b-a)*rand01();
}


// --------------------------------------------------------------------
void init( const unsigned int N )
{
    // Initialize particles
    // We will make a block of particles with a total width of 1/4 of the screen.
    float w = SIM_W/4;
    for(float y=bottom+1; y <= 10000; y+=r*.5f)
    {
        for(float x=-w; x <= w; x+=r*.5f)
        {
            if(particles.size() > N)
            {
                break;
            }

            Particle p;
            p.pos = vec2(x, y);
            p.pos_old = p.pos + 0.001f * vec2(rand01(), rand01());
            p.force = vec2(0,0);
            p.sigma = 3.f;
            p.beta = 4.f;
            particles.push_back(p);
        }
    }
}

// --------------------------------------------------------------------
void step()
{
    // UPDATE
    //
    // This modified verlet integrator has dt = 1 and calculates the velocity
    // For later use in the simulation.

    // For each particles i ...
#pragma omp parallel for
    for( int i = 0; i < (int)particles.size(); ++i )
    {
        // Apply the currently accumulated forces
        particles[i].pos += particles[i].force;

        // Restart the forces with gravity only. We'll add the rest later.
        particles[i].force = vec2( 0.0f, -::G );

        // Calculate the velocity for later.
        particles[i].vel = particles[i].pos - particles[i].pos_old;

        // If the velocity is really high, we're going to cheat and cap it.
        // This will not damp all motion. It's not physically-based at all. Just
        // a little bit of a hack.
        const float max_vel = 2.f;
        const float vel_mag = glm::length2( particles[i].vel );
        // If the velocity is greater than the max velocity, then cut it in half.
        if( vel_mag > max_vel*max_vel )
        {
            particles[i].vel *= .5f;
        }

        // Normal verlet stuff
        particles[i].pos_old = particles[i].pos;
        particles[i].pos += particles[i].vel;

        // If the Particle is outside the bounds of the world, then
        // Make a little spring force to push it back in.
        if( particles[i].pos.x < -SIM_W ) particles[i].force.x -= ( particles[i].pos.x - -SIM_W ) / 8;
        if( particles[i].pos.x >  SIM_W ) particles[i].force.x -= ( particles[i].pos.x - SIM_W ) / 8;
        if( particles[i].pos.y < bottom ) particles[i].force.y -= ( particles[i].pos.y - bottom ) / 8;
        //if( particles[i].pos.y > SIM_W * 2 ) particles[i].force.y -= ( particles[i].pos.y - SIM_W * 2 ) / 8;

        // Handle the mouse attractor.
        // It's a simple spring based attraction to where the mouse is.
        const float attr_dist2 = glm::length2( particles[i].pos - attractor );
        const float attr_l = SIM_W / 4;
        if( attracting )
        {
            if( attr_dist2 < attr_l * attr_l )
            {
                particles[i].force -= ( particles[i].pos - attractor ) / 256.0f;
            }
        }

        // Reset the nessecary items.
        particles[i].rho = 0;
        particles[i].rho_near = 0;
        particles[i].neighbors.clear();
    }

    // DENSITY
    //
    // Calculate the density by basically making a weighted sum
    // of the distances of neighboring particles within the radius of support (r)

    // For each Particle ...
#pragma omp parallel for
    for( int i = 0; i < (int)particles.size(); ++i )
    {
        particles[i].rho = 0;
        particles[i].rho_near = 0;

        // We will sum up the 'near' and 'far' densities.
        float d=0;
        float dn=0;

        // Now look at every other Particle
        for( int j = 0; j < (int)particles.size(); ++j )
        {
            // We only want to look at each pair of particles just once.
            // And do not calculate an interaction for a Particle with itself!
            if( j >= i ) continue;

            // The vector seperating the two particles
            const vec2 rij = particles[j].pos - particles[i].pos;

            // Along with the squared distance between
            const float rij_len2 = glm::length2( rij );

            // If they're within the radius of support ...
            if( rij_len2 < rsq )
            {
                // Get the actual distance from the squared distance.
                float rij_len = sqrt(rij_len2);

                // And calculated the weighted distance values
                const float q = 1 - ( rij_len / r );
                const float q2 = q * q;
                const float q3 = q2 * q;

                d += q2;
                dn += q3;

                // Accumulate on the Neighbor
                particles[j].rho += q2;
                particles[j].rho_near += q3;

                // Set up the Neighbor list for faster access later.
                Neighbor n;
                n.j = j;
                n.q = q; 
                n.q2 = q2;
                particles[i].neighbors.push_back(n);
            }
        }

        particles[i].rho += d;
        particles[i].rho_near += dn;
    }

    // PRESSURE
    //
    // Make the simple pressure calculation from the equation of state.
#pragma omp parallel for
    for( int i = 0; i < (int)particles.size(); ++i )
    {
        particles[i].press = k * ( particles[i].rho - rest_density );
        particles[i].press_near = k_near * particles[i].rho_near;
    }

    // PRESSURE FORCE
    //
    // We will force particles in or out from their neighbors
    // based on their difference from the rest density.

    // For each Particle ...
#pragma omp parallel for
    for( int i = 0; i < (int)particles.size(); ++i )
    {
        // For each of the neighbors
        vec2 dX;
        for( int ni = 0; ni < (int)particles[i].neighbors.size(); ++ni )
        {
            const Neighbor& n = particles[i].neighbors[ni];

            // The vector from Particle i to Particle j
            const vec2 rij = particles[n.j].pos - particles[i].pos;

            // calculate the force from the pressures calculated above
            const float dm 
                = n.q * ( particles[i].press + particles[n.j].press ) 
				+ n.q2 * ( particles[i].press_near + particles[n.j].press_near );

            // Get the direction of the force
            const vec2 D = glm::normalize( rij ) * dm;
            dX += D;
            particles[n.j].force += D;
        }

        particles[i].force -= dX;
    }

    // VISCOSITY
    //
    // This simulation actually may look okay if you don't compute
    // the viscosity section. The effects of numerical damping and
    // surface tension will give a smooth appearance on their own.
    // Try it.

    // For each Particle
#pragma omp parallel for
    for( int i = 0; i < (int)particles.size(); ++i )
    {
        // We'll let the color be determined by
        // ... x-velocity for the red component
        // ... y-velocity for the green-component
        // ... pressure for the blue component
        particles[i].r = 0.3f + (20 * fabs(particles[i].vel.x) );
        particles[i].g = 0.3f + (20 * fabs(particles[i].vel.y) );
        particles[i].b = 0.3f + (0.1f * particles[i].rho );

        // For each of that particles neighbors
        for( int ni = 0; ni < (int)particles[i].neighbors.size(); ++ni )
        {
            const Neighbor& n = particles[i].neighbors[ni];

            const vec2 rij = particles[n.j].pos - particles[i].pos;
            const float l = glm::length( rij );
            const float q = l / r;

            const vec2 rijn = ( rij / l );
            // Get the projection of the velocities onto the vector between them.
            const float u = glm::dot( particles[i].vel - particles[n.j].vel, rijn );
            if( u > 0 )
            {
                // Calculate the viscosity impulse between the two particles
                // based on the quadratic function of projected length.
                const vec2 I 
                    = ( 1 - q ) 
                    * ( particles[n.j].sigma * u + particles[n.j].beta * u * u ) 
                    * rijn;

                // Apply the impulses on the two particles
                particles[i].vel -= I * 0.5f;
                particles[n.j].vel += I * 0.5f;
            }
        }
    }
}

// --------------------------------------------------------------------
void display()
{
    glClearColor( 0, 0, 0, 1 );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    // create a world with dimensions x:[-SIM_W,SIM_W] and y:[0,SIM_W*2]
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho( -SIM_W, SIM_W, 0, 2*SIM_W, -1, 1 );

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // Draw Fluid Particles
    glPointSize(r*2);
    glVertexPointer( 2, GL_FLOAT, sizeof(Particle), &particles[0].pos );
    glColorPointer( 3, GL_FLOAT, sizeof(Particle), &particles[0].r );
    glEnableClientState( GL_VERTEX_ARRAY );
    glEnableClientState( GL_COLOR_ARRAY );
    glDrawArrays( GL_POINTS, 0, particles.size() );
    glDisableClientState( GL_VERTEX_ARRAY );
    glDisableClientState( GL_COLOR_ARRAY );

    glutSwapBuffers();
}

// --------------------------------------------------------------------
void idle()
{
    step();
    glutPostRedisplay();
}

// --------------------------------------------------------------------
void keyboard(unsigned char c, int x, int y)
{
    const float radius = SIM_W / 8;
    const float dtheta = 3.f;

    switch(c)
    {
        // Quit
    case 27:
    case 'q':
    case 'Q':
        exit(0);
        break;

        // If we press the space key, add some particles.
    case ' ':
        for( float y = SIM_W * 2 - radius; y <= SIM_W * 2 + radius; y += r * .5f )
        {
            for( float x = -radius; x <= radius; x += r * .5f )
            {
                Particle p;
                p.pos = p.pos_old = vec2(x , y) + vec2(rand01(), rand01());
                p.force = vec2(0,0);
                p.sigma = 3.f;
                p.beta = 4.f;

                if( glm::length2( p.pos - vec2( 0, SIM_W * 2 ) ) < radius * radius )
                {
                    particles.push_back(p);
                }
            }
        }
        break;
    }
}

// --------------------------------------------------------------------
void motion(int x, int y)
{
    // This simply updates the location of the mouse attractor.
    int window_w = glutGet( GLUT_WINDOW_WIDTH );
    int window_h = glutGet( GLUT_WINDOW_HEIGHT );
    float relx = (float)(x - window_w/2) / window_w;
    float rely = -(float)(y - window_h) / window_h;
    vec2 mouse = vec2(relx*SIM_W*2, rely*SIM_W*2);
    attractor = mouse;
}

// --------------------------------------------------------------------
void mouse(int button, int state, int x, int y)
{
    if(state == GLUT_DOWN)
    {
        attracting = true;
    }
    else
    {
        attracting = false;
        attractor = vec2(SIM_W * 99, SIM_W * 99);
    }
}

// --------------------------------------------------------------------
int main( int argc, char** argv )
{
#if 1
    ofstream file( "benchmark.txt" );
    typedef iostreams::tee_device< std::ostream, std::ofstream > Tee;
    typedef iostreams::stream< Tee > TeeStream;
    Tee tee( std::cout, file );
    TeeStream log( tee );

    const int steps = 3000;
    log << "--------------------------------" << endl;
    log << "Number of steps: " << steps << endl;
    for( unsigned int size = 10; size <= 13; ++size )
    {
        const unsigned int count = ( 1 << size );
        log << "Number of particles: " << count << endl;

        init( count );

        typedef chrono::high_resolution_clock clock_t;
        clock_t::time_point beg, end;

        beg = clock_t::now();
        for( int i = 0; i < steps; ++i )
        {
            step();
        }
        end = clock_t::now();

        log << "Elapsed time: " << chrono::duration_cast< chrono::milliseconds >( end - beg ) << endl;
        log << "Microseconds per step: " << chrono::duration_cast< chrono::microseconds >( end - beg ).count() / (double)steps << endl;
        log << endl;
    }

    return 0;
#endif

    init( 2048 );

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE);
    glutInitWindowSize(512, 512);
    glutCreateWindow("SPH");

    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutIdleFunc(idle);
    glutMotionFunc(motion);
    glutMouseFunc(mouse);

    glutMainLoop();
}