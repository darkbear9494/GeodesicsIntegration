#include <iostream>
#include <omp.h> // OpenMP编程需要包含的头文件

int main()
{
    std::cout << omp_get_num_procs() << std::endl;

#pragma omp parallel //指明下面大括号内部为并行区域
    {
        std::cout << "OpenMP" << std::endl;
    }

    return 0;
}
