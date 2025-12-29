#ifndef __SIEVE2_C__
#define __SIEVE2_C__

#include "include.h"

// Follow instructions based on HW3 guide. Thanks, TA.
void sieve2(unsigned long long *global_count, unsigned long long n, int pnum, int pid)
{
    unsigned long long low_value = 3 + 2 * ((pid * ((n - 3) / 2 + 1)) / pnum); // smallest value handled by this process
    unsigned long long high_value = 3 + 2 * (((pid + 1) * ((n - 3) / 2 + 1)) / pnum) - 2; // largest value handled by this process
    unsigned long long size = (high_value - low_value) / 2 + 1; // number of integers handled by this process

    if (1 + (n - 1) / pnum < (int)sqrt((double)n)) // high_value of process 0 should be larger than floor(sqrt(n))
    {
        if (pid == 0)
            printf("Error: Too many processes.\n");
        MPI_Finalize();
        exit(0);
    }

    char *marked = (char*)malloc(size); // array for marking multiples. 1 means multiple and 0 means prime
    if (marked == NULL)
    {
        printf("Error: Cannot allocate enough memory.\n");
        MPI_Finalize();
        exit(0);
    }
    memset(marked, 0, size);

    unsigned long long sqrt_n = (unsigned long long)sqrt((double)n);
    unsigned long long sieve_size = (sqrt_n - 3) / 2 + 1; // size for local sieving primes
    char *sieve = (char*)malloc(sieve_size); // array for local sieving primes
    if (sieve == NULL)
    {
        printf("Error: Cannot allocate enough memory for sieve.\n");
        MPI_Finalize();
        exit(0);
    }
    memset(sieve, 0, sieve_size);

    // Generate sieving primes locally
    for (unsigned long long i = 0; i < sieve_size; i++)
    {
        if (sieve[i] == 0) // Found a prime
        {
            unsigned long long prime = 2 * i + 3;
            for (unsigned long long j = (prime * prime - 3) / 2; j < sieve_size; j += prime)
                sieve[j] = 1;
        }
    }

    // Use sieving primes to mark composites in the range
    for (unsigned long long i = 0; i < sieve_size; i++)
    {
        if (sieve[i] == 0) // Found a sieving prime
        {
            unsigned long long prime = 2 * i + 3;
            unsigned long long first; // First multiple to mark
            if (prime * prime > low_value)
                first = ((prime * prime) - low_value) / 2;
            else if (low_value % prime == 0)
                first = 0;
            else
                first = ((prime * ((((low_value / prime) + 1) / 2) * 2 + 1)) - low_value) / 2;

            for (unsigned long long j = first; j < size; j += prime)
                marked[j] = 1;
        }
    }

    unsigned long long count = 0; // Local count of primes
    for (unsigned long long i = 0; i < size; i++)
        if (marked[i] == 0)
            count++;

    if (pid == 0) count++; // Add prime 2 to the count
    MPI_Reduce(&count, global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Free allocated memory
    free(marked);
    free(sieve);
}

#endif
