#ifndef __SIEVE3_C__
#define __SIEVE3_C__

#include "include.h"

// Follow instructions based on HW3 guide and lecture slides.
#define BLOCK_SIZE 1024 * 32  // Based on cache optimization hints

void sieve3(unsigned long long *global_count, unsigned long long n, int pnum, int pid)
{
    unsigned long long sqrt_n = (unsigned long long)sqrt((double)n);

    // Step 1: Calculate local primes up to sqrt(n)
    unsigned long long local_sieve_size = (sqrt_n - 3) / 2 + 1; // Only odd numbers up to sqrt(n)
    char *local_sieve = (char *)malloc(local_sieve_size);
    if (!local_sieve)
    {
        printf("Error: Unable to allocate memory for local sieve array.\n");
        MPI_Finalize();
        exit(0);
    }
    memset(local_sieve, 0, local_sieve_size);

    // Generate local primes
    for (unsigned long long i = 0; i < local_sieve_size; i++)
    {
        if (local_sieve[i] == 0)
        {
            unsigned long long prime = 2 * i + 3;
            for (unsigned long long j = (prime * prime - 3) / 2; j < local_sieve_size; j += prime)
                local_sieve[j] = 1;
        }
    }

    // Store local primes for marking later
    unsigned long long *local_primes = (unsigned long long *)malloc(local_sieve_size * sizeof(unsigned long long));
    unsigned long long local_prime_count = 0;
    for (unsigned long long i = 0; i < local_sieve_size; i++)
    {
        if (local_sieve[i] == 0)
            local_primes[local_prime_count++] = 2 * i + 3;
    }
    free(local_sieve);

    // Step 2: Define the range for this process
    unsigned long long low_value = 3 + 2 * ((pid * ((n - 3) / 2 + 1)) / pnum);
    unsigned long long high_value = 3 + 2 * (((pid + 1) * ((n - 3) / 2 + 1)) / pnum) - 2;
    unsigned long long size = (high_value - low_value) / 2 + 1;

    // Allocate memory for marking composites
    char *marked = (char *)malloc(BLOCK_SIZE);
    if (!marked)
    {
        printf("Error: Unable to allocate memory for marked array.\n");
        MPI_Finalize();
        exit(0);
    }

    unsigned long long count = 0;

    // Step 3: Block processing for cache efficiency
    for (unsigned long long block_low = low_value; block_low <= high_value; block_low += 2 * BLOCK_SIZE)
    {
        unsigned long long block_high = block_low + 2 * BLOCK_SIZE - 1;
        if (block_high > high_value)
            block_high = high_value;

        unsigned long long block_size = (block_high - block_low) / 2 + 1;
        memset(marked, 0, block_size);

        // Mark multiples of local primes in the current block
        for (unsigned long long i = 0; i < local_prime_count; i++)
        {
            unsigned long long prime = local_primes[i];
            unsigned long long first;

            if (prime * prime > block_low)
                first = (prime * prime - block_low) / 2;
            else if (block_low % prime == 0)
                first = 0;
            else
            {
                unsigned long long temp = (block_low + prime) / prime;
                if (temp % 2 == 0)
                    temp++;
                first = (temp * prime - block_low) / 2;
            }

            for (unsigned long long j = first; j < block_size; j += prime)
                marked[j] = 1;
        }

        // Count primes in the current block
        for (unsigned long long i = 0; i < block_size; i++)
            if (marked[i] == 0)
                count++;
    }

    // Include prime 2 in the count
    if (pid == 0) count++;

    // Step 4: Reduce the counts across all processes
    MPI_Reduce(&count, global_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    free(marked);
    free(local_primes);
}

#endif
