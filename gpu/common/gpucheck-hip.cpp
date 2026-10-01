// gpucheck-hip <expected GPUs>: HIP counterpart of gpucheck.cu for Frontier
// (submit-bundle-frontier.sh). Creates a context on every GPU of the node and
// prints one line, "<host> OK <n> GPUs" or "<host> BAD ...", so a node whose
// GPUs are missing or fail context creation is left out instead of killing
// the job. Always exits 0 so srun keeps the other tasks.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

#include <hip/hip_runtime.h>

int main(int argc, char** argv)
{
    const int expected = argc > 1 ? atoi(argv[1]) : 1;
    char host[256];
    gethostname(host, sizeof(host));
    if (char* dot = strchr(host, '.')) {
        *dot = 0;
    }

    int n = 0;
    hipError_t err = hipGetDeviceCount(&n);
    if (err != hipSuccess) {
        printf("%s BAD hipGetDeviceCount: %s\n", host, hipGetErrorString(err));
        return 0;
    }
    if (n < expected) {
        printf("%s BAD %d of %d GPUs\n", host, n, expected);
        return 0;
    }
    for (int d = 0; d < n; d++) {
        err = hipSetDevice(d);
        if (err == hipSuccess) {
            err = hipFree(0);
        }
        if (err != hipSuccess) {
            printf("%s BAD GPU %d: %s\n", host, d, hipGetErrorString(err));
            return 0;
        }
    }
    printf("%s OK %d GPUs\n", host, n);
    return 0;
}
