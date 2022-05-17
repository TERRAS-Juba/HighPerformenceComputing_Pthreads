#include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>
#include <pthread.h>
// --------------------------------------------
// 1)- Declaration des structures
// --------------------------------------------
#define NUM_THREADS 8
pthread_mutex_t mutex_test;
int nd=3;
int np=50;
int step_num=10;
double dt=0.1;
int seed = 123456789;
double mass = 1.0;
double kinetic;
double potential;
typedef struct thread_strcut {
    int id;
    int dimension;
    int nb_iterations;
    int debut;
    double *pos;
    double *vel;
    double *acc;
    double *force;
};
// --------------------------------------------
// 2)- Declaration des entêtes des fonctions
// --------------------------------------------
int main ( int argc, char *argv[] );
void compute ( int np, int nd, double pos[], double vel[],
               double mass, double f[], double *pot, double *kin );
double cpu_time ( );
double dist ( int nd, double r1[], double r2[], double dr[] );
void initialize ( int np, int nd, double pos[], double vel[], double acc[] );
void *initialize_parallele ( void *t );
void r8mat_uniform_ab ( int m, int n, double a, double b, int *seed, double r[] );
void timestamp ( );
void update ( int np, int nd, double pos[], double vel[], double f[],
              double acc[], double mass, double dt );
// --------------------------------------------
// 3)- Programme principale
// --------------------------------------------
int main ( int argc, char *argv[] )
{
    double ctime;
    double e0;
    // Matrice pour le programme sequentielle
    double *acc;
    double *vel;
    double *force;
    double *pos;
    // Matrice pour le programme parallele
    double *acc_parallele;
    double *vel_parallele;
    double *force_parallele;
    double *pos_parallele;
    // Variables pour les iterations
    int step;
    int step_print;
    int step_print_index;
    int step_print_num;
    timestamp ( );
    printf ( "\n" );
    printf ( "MD\n" );
    printf ( "  C version\n" );
    printf ( "  A molecular dynamics program.\n" );
    // -----------------------------------------------------------
    // 3-1)- Introduction des parameters via la ligne de commande
    // -----------------------------------------------------------
    /*if ( 1 < argc )
    {
        nd = atoi ( argv[1] );
    }
    else
    {
        printf ( "\n" );
        printf ( "  Enter ND, the spatial dimension (2 or 3).\n" );
        scanf ( "%d", &nd );
    }
    if ( 2 < argc )
    {
        np = atoi ( argv[2] );
    }
    else
    {
        printf ( "\n" );
        printf ( "  Enter NP, the number of particles (500, for instance).\n" );
        scanf ( "%d", &np );
    }
    if ( 3 < argc )
    {
        step_num = atoi ( argv[3] );
    }
    else
    {
        printf ( "\n" );
        printf ( "  Enter ND, the number of time steps (500 or 1000, for instance).\n" );
        scanf ( "%d", &step_num );
    }
    if ( 4 < argc )
    {
        dt = atof ( argv[4] );
    }
    else
    {
        printf ( "\n" );
        printf ( "  Enter DT, the size of the time step (0.1, for instance).\n" );
        scanf ( "%lf", &dt );
    }
    printf ( "\n" );
    printf ( "  ND, the spatial dimension, is %d\n", nd );
    printf ( "  NP, the number of particles in the simulation, is %d\n", np );
    printf ( "  STEP_NUM, the number of time steps, is %d\n", step_num );
    printf ( "  DT, the size of each time step, is %f\n", dt );*/
    // -----------------------------------------------------------
    // 3-2)- Allocation memoire
    // -----------------------------------------------------------
    acc = ( double * ) malloc ( nd * np * sizeof ( double ) );
    force = ( double * ) malloc ( nd * np * sizeof ( double ) );
    pos = ( double * ) malloc ( nd * np * sizeof ( double ) );
    vel = ( double * ) malloc ( nd * np * sizeof ( double ) );
    pos_parallele = ( double * ) malloc ( nd * np * sizeof ( double ) );
    force_parallele = ( double * ) malloc ( nd * np * sizeof ( double ) );
    acc_parallele = ( double * ) malloc ( nd * np * sizeof ( double ) );
    vel_parallele = ( double * ) malloc ( nd * np * sizeof ( double ) );
    // -----------------------------------------------------------
    // 3-3)- Affichage
    // -----------------------------------------------------------
    printf ( "\n" );
    printf ( "  At each step, we report the potential and kinetic energies.\n" );
    printf ( "  The sum of these energies should be a constant.\n" );
    printf ( "  As an accuracy check, we also print the relative error\n" );
    printf ( "  in the total energy.\n" );
    printf ( "\n" );
    printf ( "      Step      Potential       Kinetic        (P+K-E0)/E0\n" );
    printf ( "                Energy P        Energy K       Relative Energy Error\n" );
    printf ( "\n" );
    step_print = 0;
    step_print_index = 0;
    step_print_num = 10;
    // -----------------------------------------------------------
    // 3-4)- Portion parallele
    // -----------------------------------------------------------
    ctime = cpu_time ( ); // Time CPU
    // Les structures
    pthread_t threads[NUM_THREADS];
    pthread_mutex_init(&mutex_test,NULL);
    struct thread_strcut threads_structs[NUM_THREADS];
    // Remplissage des structures de threads
    for (int i = 0; i < NUM_THREADS ; ++i) {
        threads_structs[i].id=i;
        threads_structs[i].dimension=nd;
        threads_structs[i].nb_iterations=np/NUM_THREADS;
        threads_structs[i].debut=i*threads_structs[i].nb_iterations;
        threads_structs[i].pos=pos_parallele;
        threads_structs[i].vel=vel_parallele;
        threads_structs[i].acc=acc_parallele;
        threads_structs[i].force=force_parallele;
    }
    // Creation des threads
    for (int j = 0; j <= step_num; ++j) {
        if (j == 0) {
            for (int i = 0; i < NUM_THREADS; ++i) {
                pthread_create(&threads[i], NULL, initialize_parallele, (void *) &threads_structs[i]);
                for (i = 0; i < NUM_THREADS; ++i) {
                    pthread_join(threads[i], NULL);
                }
            }
        }
    }
    // -----------------------------------------------------------
    // 3-4)- Portion sequentielle
    // -----------------------------------------------------------
    for ( step = 0; step <= step_num; step++ )
    {
        if ( step == 0 )
        {
            initialize ( np, nd, pos, vel, acc );
        }
        else
        {
            update ( np, nd, pos, vel, force, acc, mass, dt );
        }

        compute ( np, nd, pos, vel, mass, force, &potential, &kinetic );

        if ( step == 0 )
        {
            e0 = potential + kinetic;
        }

        if ( step == step_print )
        {
            printf ( "  %8d  %14f  %14f  %14e\n", step, potential, kinetic,
                     ( potential + kinetic - e0 ) / e0 );
            step_print_index = step_print_index + 1;
            step_print = ( step_print_index * step_num ) / step_print_num;
        }

    }
    // -----------------------------------------------------------
    // 3-5)- Affichage du temps d'execution seuqentielle
    // -----------------------------------------------------------
    ctime = cpu_time ( ) - ctime;
    printf ( "\n" );
    printf ( "  Elapsed cpu time: %f seconds.\n", ctime );
    for (int i = 0; i < np; ++i) {
        for (int j = 0; j < nd; ++j) {
            printf("%f |",pos_parallele[j+i*nd]);
        }
    }
    printf("\n-----------------------------------------------------------\n");
    for (int i = 0; i < np; ++i) {
        for (int j = 0; j < nd; ++j) {
            printf("%f |",pos[j+i*nd]);
        }
    }
    // -----------------------------------------------------------
    // 3-6)- Liberation des ressources en memoire
    // -----------------------------------------------------------
    free ( acc );
    free ( force );
    free ( pos );
    free ( vel );
    // -----------------------------------------------------------
    // 3-7)- Terminaison du programme
    // -----------------------------------------------------------
    printf ( "\n" );
    printf ( "MD\n" );
    printf ( "  Normal end of execution.\n" );
    printf ( "\n" );
    timestamp ( );
    return 0;
}
/******************************************************************************/
void compute ( int np, int nd, double pos[], double vel[], double mass,
               double f[], double *pot, double *kin )
{
    double d;
    double d2;
    int i;
    int j;
    int k;
    double ke;
    double pe;
    double PI2 = 3.141592653589793 / 2.0;
    double rij[3];

    pe = 0.0;
    ke = 0.0;

    for ( k = 0; k < np; k++ )
    {
/*
  Compute the potential energy and forces.
*/
        for ( i = 0; i < nd; i++ )
        {
            f[i+k*nd] = 0.0;
        }

        for ( j = 0; j < np; j++ )
        {
            if ( k != j )
            {
                d = dist ( nd, pos+k*nd, pos+j*nd, rij );
/*
  Attribute half of the potential energy to particle J.
*/
                if ( d < PI2 )
                {
                    d2 = d;
                }
                else
                {
                    d2 = PI2;
                }

                pe = pe + 0.5 * pow ( sin ( d2 ), 2 );

                for ( i = 0; i < nd; i++ )
                {
                    f[i+k*nd] = f[i+k*nd] - rij[i] * sin ( 2.0 * d2 ) / d;
                }
            }
        }
/*
  Compute the kinetic energy.
*/
        for ( i = 0; i < nd; i++ )
        {
            ke = ke + vel[i+k*nd] * vel[i+k*nd];
        }
    }

    ke = ke * 0.5 * mass;

    *pot = pe;
    *kin = ke;

    return;
}
/*******************************************************************************/
double cpu_time ( )
{
    double value;

    value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

    return value;
}
/******************************************************************************/
double dist ( int nd, double r1[], double r2[], double dr[] )
{
    double d;
    int i;

    d = 0.0;
    for ( i = 0; i < nd; i++ )
    {
        dr[i] = r1[i] - r2[i];
        d = d + dr[i] * dr[i];
    }
    d = sqrt ( d );

    return d;
}
/******************************************************************************/
void initialize ( int np, int nd, double pos[], double vel[], double acc[] )
{
    int i;
    int j;
    int seed;
/*
  Set positions.
*/
    seed = 123456789;
    r8mat_uniform_ab ( nd, np, 0.0, 10.0, &seed, pos );
/*
  Set velocities.
*/
    for ( j = 0; j < np; j++ )
    {
        for ( i = 0; i < nd; i++ )
        {
            vel[i+j*nd] = 0.0;
        }
    }
/*
  Set accelerations.
*/
    for ( j = 0; j < np; j++ )
    {
        for ( i = 0; i < nd; i++ )
        {
            acc[i+j*nd] = 0.0;
        }
    }

    return;
}
/******************************************************************************/
void *initialize_parallele ( void *t )
{
    struct thread_strcut *structure=(struct thread_strcut*)t;
    int i;
    int j;
/*
  Set positions.
*/
    if(structure->id==0){
        r8mat_uniform_ab ( nd, np, 0.0, 10.0, &seed, structure->pos );
    }
/*
  Set velocities.
*/
    for ( j = structure->debut; j < structure->debut+structure->nb_iterations; j++ )
    {
        for ( i = 0; i < structure->dimension; i++ )
        {
            structure->vel[i+j*structure->dimension] = 0.0;
        }
    }
/*
  Set accelerations.
*/
    for ( j = structure->debut; j < structure->debut+structure->nb_iterations; j++ )
    {
        for ( i = 0; i < structure->dimension; i++ )
        {
            structure->acc[i+j*structure->dimension] = 0.0;
        }
    }
    pthread_exit(NULL);
}
/******************************************************************************/
void r8mat_uniform_ab ( int m, int n, double a, double b, int *seed, double r[] )
{
    int i;
    const int i4_huge = 2147483647;
    int j;
    int k;

    if ( *seed == 0 )
    {
        fprintf ( stderr, "\n" );
        fprintf ( stderr, "R8MAT_UNIFORM_AB - Fatal error!\n" );
        fprintf ( stderr, "  Input value of SEED = 0.\n" );
        exit ( 1 );
    }

    for ( j = 0; j < n; j++ )
    {
        for ( i = 0; i < m; i++ )
        {
            k = *seed / 127773;

            *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

            if ( *seed < 0 )
            {
                *seed = *seed + i4_huge;
            }
            r[i+j*m] = a + ( b - a ) * ( double ) ( *seed ) * 4.656612875E-10;
        }
    }

    return;
}
/******************************************************************************/
void timestamp ( )
{
# define TIME_SIZE 40

    static char time_buffer[TIME_SIZE];
    const struct tm *tm;
    time_t now;

    now = time ( NULL );
    tm = localtime ( &now );

    strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

    printf ( "%s\n", time_buffer );

    return;
# undef TIME_SIZE
}
/******************************************************************************/

void update ( int np, int nd, double pos[], double vel[], double f[],
              double acc[], double mass, double dt )

{
    int i;
    int j;
    double rmass;

    rmass = 1.0 / mass;

    for ( j = 0; j < np; j++ )
    {
        for ( i = 0; i < nd; i++ )
        {
            pos[i+j*nd] = pos[i+j*nd] + vel[i+j*nd] * dt + 0.5 * acc[i+j*nd] * dt * dt;
            vel[i+j*nd] = vel[i+j*nd] + 0.5 * dt * ( f[i+j*nd] * rmass + acc[i+j*nd] );
            acc[i+j*nd] = f[i+j*nd] * rmass;
        }
    }

    return;
}