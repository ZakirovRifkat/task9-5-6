#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
void gauss(double** mas1, double** mas2, int n)
{
	int k = 0, ind;
	while (k < n)
	{
		double t = mas1[k][k];
		if (t == 0)
		{
			for (int i = 0; i < n; i++)
				if (mas1[i][k] != 0)
				{
					ind = i;
					break;
				}
			for (int i = 0; i < n; i++)
			{
				double t = mas1[k][i], t1 = mas2[k][i];
				mas1[k][i] = mas1[ind][i], mas2[k][i] = mas2[ind][i];
				mas1[ind][i] = t, mas2[ind][i] = t1;
			}
			t = mas1[k][k];
		}
		for (int i = 0; i < n; i++)
		{
			mas1[k][i] /= t;
			mas2[k][i] /= t;
		}
		for (int i = k + 1; i < n; i++)
		{
			double t = mas1[i][k];
			for (int j = 0; j < n; j++)
			{
				mas1[i][j] -= mas1[k][j] * t;
				mas2[i][j] -= mas2[k][j] * t;
			}
		}
		k++;
	}
	for (k = n - 1; k > 0; k--)
	{
		for (int j = k - 1; j >= 0; j--)
		{
			double t = mas1[j][k];
			mas1[j][k] -= mas1[k][k] * t;
			for (int i = n - 1; i >= 0; i--)
			{
				mas2[j][i] -= mas2[k][i] * t;
			}
		}
	}
}



int main()
{ int i,j,k,n;
    float A[100][100],c,x[10],sum=0.0;
	double** mas1, ** mas2, ** mas3;
    setlocale(LC_ALL,"Russian");
printf("1-Решить систему линейных уравнений методом Гаусса(n уравнений, n неизвестных)\n");
printf("2-Найти обратную матрицу методом Гаусса (сделать проверку)\n");
printf("|Для выхода введите 0|\n\n");
int p=-1;
while (p!=0)
{
    printf("Выберите задание:");
    scanf("%d",&p);
    switch(p){
case 1:

    printf("\nВведите кол-во переменных ");
    scanf("%d",&n);
    printf("\nВведите элементы матрицы:\n\n");
    for(i=1; i<=n; i++)
    {
        for(j=1; j<=(n+1); j++)
        {
            printf("A[%d][%d] : ", i,j);
            scanf("%f",&A[i][j]);
        }
    }
    for(j=1; j<=n; j++)
    {
        for(i=1; i<=n; i++)
        {
            if(i>j)
            {
                c=A[i][j]/A[j][j];
                for(k=1; k<=n+1; k++)
                {
                    A[i][k]=A[i][k]-c*A[j][k];
                }
            }
        }
    }
    x[n]=A[n][n+1]/A[n][n];

    for(i=n-1; i>=1; i--)
    {
        sum=0;
        for(j=i+1; j<=n; j++)
        {
            sum=sum+A[i][j]*x[j];
        }
        x[i]=(A[i][n+1]-sum)/A[i][i];
    }
    printf("\nОтвет: \n");
    for(i=1; i<=n; i++)
    {
        printf("\nx%d=%f\t",i,x[i]);
    }
    printf("\n\n");
    break;

    case 2:
        printf("Введите размер матрицы: ");
	scanf("%d", &n);

	mas1 = (double**)malloc(n * sizeof(double*));
	mas2 = (double**)malloc(n * sizeof(double*));
	mas3 = (double**)malloc(n * sizeof(double*));


	for (int i = 0; i < n; ++i) {
		mas1[i] = (double*)malloc(n * sizeof(double));
		mas2[i] = (double*)malloc(n * sizeof(double));
		mas3[i] = (double*)malloc(n * sizeof(double));
	}
	printf("Элементы матрицы: \n");

	for (int i = 0; i < n; i++)

		for (int j = 0; j < n; j++)

		{
			if (i == j)
				mas2[i][j] = 1;
			else mas2[i][j] = 0;

			scanf("%lf", &mas1[i][j]);
			mas3[i][j] = mas1[i][j];
		}

	gauss(mas1, mas2, n);
	printf("\nОтвет: \n");
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			printf("%lf ", mas2[i][j]);
		printf("\n");
	}

	double** mas4;
	mas4 = (double**)malloc(n * sizeof(double*));
	for (int i = 0; i < n; i++)
	{
		mas4[i] = (double*)malloc(n * sizeof(double));

		for (int j = 0; j < n; j++)

		{
			mas4[i][j] = 0;

			for (int k = 0; k < n; k++)
				mas4[i][j] += mas3[i][k] * mas2[k][j];

			if (fabs(mas4[i][j]) < 0.000001)
				mas4[i][j] = 0;
		}

	}
	printf("\n");

    break;
    }


}
    return 0;
}
