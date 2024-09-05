#include <iostream>
using namespace std;
#include <cmath>
#include <cstdlib>

// I just converted my python code for N-body to c++
// And then modified it. 
double m_ast =10*1.989e18;
double m1 = 1.989e30;
double m2 = 1.898e27;
double rj = 7.78479e11;
double G = 6.674e-11;
double pi = 3.1415926;
double w = sqrt(G*(m1+m2)/pow(rj,3));
double r1 = -m2*rj/(m1+m2); 
double r2 = m1*rj/(m1+m2);
struct Particle
{
    double mass;
    double x, y;
    double vx, vy;

};

double drand() {
    return (double)rand()/(double)RAND_MAX;
}

// Generate x,y,vx,vy for asteroids
Particle gen()
{
    Particle a;
    a.mass = m_ast;
    double r = 2.7e11 + 2.7e11 * drand();
    double th = pi;
    th *= 2*drand();
    double root = m1;
    root *= G;
    root /= r;
    double v = sqrt(root) - w * r;
    v += v * (2 * drand() - 1) * 0.05;
    double vth = th * (1 + (2 * drand() - 1) * 0.05);
    a.x = r * cos(th);
    a.y = r * sin(th);
    a.vx = -v * sin(vth);
    a.vy = v * cos(vth);
    return a;
}

void update(Particle particles[], int numParticles, double dt)
{
    for (int i = 0; i < 2; i++) // To neglect the asteroid-asteroid interaction, 'i' will only take value 0 and 1. 
    {
        Particle &particle1 = particles[i];
        for (int j = i + 1; j < numParticles; j++)
        {
            Particle &particle2 = particles[j];

            double dx = particle2.x - particle1.x;
            double dy = particle2.y - particle1.y;
            double distance = sqrt(dx * dx + dy * dy);

            double f_x = (G * particle1.mass * particle2.mass * dx) / (distance * distance * distance);
            double f_y = (G * particle1.mass * particle2.mass * dy) / (distance * distance * distance);

            particle1.vx += f_x / particle1.mass * dt;
            particle1.vy += f_y / particle1.mass * dt;

            particle2.vx -= f_x / particle2.mass * dt;
            particle2.vy -= f_y / particle2.mass * dt;
        }
    }
    for (int i = 0; i < numParticles; i++)
    {
        Particle &particle = particles[i];
        particle.x += particle.vx * dt;
        particle.y += particle.vy * dt;
    }
}

void simulate(Particle particles[], int numParticles, double dt, int steps)
{
    double *x[numParticles];
    double *y[numParticles];
    for (int i = 0; i < numParticles; i++)
    {
        x[i] = new double[steps];
        y[i] = new double[steps];
    }

    for (int step = 0; step < steps; step++)
    {
        for (int i = 0; i < numParticles; i++)
        {
            x[i][step] = particles[i].x;
            y[i][step] = particles[i].y;
        }
        update(particles, numParticles, dt);
    }

    for (int i = 0; i < numParticles; i++)
    {

        for (int step = 0; step < steps; step++)
        {
            cout << "(" << x[i][step] << ", " << y[i][step] << ") ";
        }
        cout << endl;

        // Clean up memory
        delete[] x[i];
        delete[] y[i];
    }
}

int main()
{
    Particle particles[10000];
    particles[0] = Particle();
    particles[1] = Particle();
    particles[0].mass = 1.989e30;
    particles[0].x = 0;
    particles[0].y = 0;
    particles[0].vx = 0;
    particles[0].vy = 0;
    particles[1].mass = 1.8982e27;
    particles[1].x = 5.2038 * 1.496e11;
    particles[1].y = 0;
    particles[1].vx = 0;
    particles[1].vy = 1.307e4;
    for (int i = 2; i < 10000; i++)
    {
        particles[i] = gen();
    }

    double dt = 3600;  // time step in seconds
    int steps = 43800000; // number of simulation steps

    simulate(particles, sizeof(particles) / sizeof(particles[0]), dt, steps);

    return 0;
}

