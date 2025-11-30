#include "onyx/app/app.hpp"
#include "tkit/multiprocessing/thread_pool.hpp"
#include "tkit/profiling/macros.hpp"

#define GENE_MAX_WORKERS (ONYX_MAX_THREADS - 1)

void RunApp()
{
    Onyx::Window::Specs spc;
    spc.Name = "Physics Eugene";

    Onyx::SingleWindowApp app{spc};
    app.Run();
}

int main(int argc, char **argv)
{
    TKIT_PROFILE_NOOP();
    TKit::ThreadPool threadPool{GENE_MAX_WORKERS};
    Onyx::Core::Initialize(Onyx::Specs{.TaskManager = &threadPool});
    RunApp();
    Onyx::Core::Terminate();
}
