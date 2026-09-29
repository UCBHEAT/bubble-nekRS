// gpucheck <expected GPUs>: creates a CUDA context on every GPU of the node
// and prints one line, "<host> OK <n> GPUs" or "<host> BAD ...".
// submit-bundle.sh runs it once per node so a node whose GPUs are missing or
// fail context creation (CUDA error 999) is left out instead of killing the
// job. Always exits 0 so mpiexec keeps the other ranks.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

#include <cuda_runtime.h>

int main(int argc, char** argv)
{
    const int expected = argc > 1 ? atoi(argv[1]) : 1;
    char host[256];
    gethostname(host, sizeof(host));
    if (char* dot = strchr(host, '.')) {
        *dot = 0;
    }

    int n = 0;
    cudaError_t err = cudaGetDeviceCount(&n);
    if (err != cudaSuccess) {
        printf("%s BAD cudaGetDeviceCount: %s\n", host, cudaGetErrorString(err));
        return 0;
    }
    if (n < expected) {
        printf("%s BAD %d of %d GPUs\n", host, n, expected);
        return 0;
    }
    for (int d = 0; d < n; d++) {
        err = cudaSetDevice(d);
        if (err == cudaSuccess) {
            err = cudaFree(0);
        }
        if (err != cudaSuccess) {
            printf("%s BAD GPU %d: %s\n", host, d, cudaGetErrorString(err));
            return 0;
        }
    }
    printf("%s OK %d GPUs\n", host, n);
    return 0;
}
