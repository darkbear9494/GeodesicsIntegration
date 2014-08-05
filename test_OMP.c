#include <stdio.h>
#include <omp.h>

int main(){
	int i;
	printf("%d\n", omp_get_num_procs());
#pragma omp parallel
	{
	printf("ok\n");
#pragma omp for
	for(i = 0; i < 10; i++){
		printf("%d\n", i);
	}

	}
}
/*
int main(){
	int i;
	#pragma omp parallel for
	for(i = 0; i < 10; i++){
		printf("%d\n", i);
	}
	return;
}
*/
