#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>



int rand(void);
float aleatorio(void);
float powf(float x, float y);
int poblar(int *red, float p, int dim);
int imprimir(int *red, int dim);
int masa(int *red, int dim,int *etiqueta_percolante, int *mass);
int problema2a(int dim, int N);
int problema2b(int N);
int campo_med(int *red,int dim,int i);
float M(int *red, int dim);
int energia_con_int(int *red, int dim, float B,float J);
float energia_sin_int(int *red, int dim, float B);
int p_libre(float B,float *p_l);
int p_int(float J,float *p_i);
int flip_spin(int *red,int dim,int *spin,float *delta_E,float J,float B,float *p_l,float *p_i);
int flip_libre(int *red,int dim,int *spin,float *delta_E,float B,float *p_l,int *indice);


int main(int argc,char *argv[]){
	int iteraciones,dim;


	double total_time;
	clock_t start, end;
	start = clock();

	sscanf(argv[1], "%d", & dim);
	sscanf(argv[2], "%d", & iteraciones);
	
	srand(time(NULL)); 



	problema2a(dim,iteraciones);

	end = clock();								//time count stops 


	total_time = ((double) (((float)(end - start)) / (float)CLOCKS_PER_SEC)/60.0);	
	printf("\nEl tiempo (minutos) requerido es:  %f \n", total_time);		//calulate total time

	
	return 0;
}





float aleatorio(void){
	float a;
	a=((float)rand())/((float)RAND_MAX);
return a;
}



int p_int(float J,float *p_i){
	int i;
	for(i=0;i<5;i++){
		*(p_i+i)=exp(-4*(1.0-0.5*((float)(i)))*J);
	}
return 0;
}

int p_libre(float B,float *p_l){
	int i;
	for(i=0;i<2;i++){
		*(p_l+i)=exp(2.0*(1.0-2.0*((float)i))*B);
	}
return 0;
}

int campo_med(int *red,int dim,int i){
	int cmed,j;
	
	if(i>dim && i<dim*dim-dim && i%dim!=0 && i%(dim-1)!=0){                          //Los casos en los que no estan en el borde de la red 
		cmed=*(red+i+1)+*(red+i-1)+*(red+i+dim)+*(red+i-dim);
	}
	else{
		if(i<dim){								//borde superior
			if(i==(dim-1) || i==0){
				if(i==(dim-1)){						//laterales del borde superior
					cmed=*(red)+*(red+i-1)+*(red+i+dim)+*(red+i+dim*dim-dim);}	//lateral izq
				else{cmed=*(red+1)+*(red+dim-1)+*(red+i+dim)+*(red+dim*dim-dim);}}	//lateral der
			else{
				cmed=*(red+i+1)+*(red+i-1)+*(red+i+dim)+*(red+i+dim*dim-dim);	//borde superior sin laterales
			}
		}
		else{										//borde inferior
			if(i==(dim*dim-1) || i==(dim*dim-dim)){
				if(i==(dim*dim-1)){						//laterales del borde inferior
					cmed=*(red+dim*dim-dim)+*(red+i-1)+*(red+i-(dim*dim-dim))+*(red+i-dim);}     //lateral izq
				else{cmed=*(red+i+1)+*(red+dim*dim-1)+*(red+i-(dim*dim-dim))+*(red+dim*dim-dim);}}	//laterla der
			else{
				cmed=*(red+i+1)+*(red+i-1)+*(red+i-(dim*dim-dim))+*(red+i-dim);	//borde inferior sin laterales
			}
			
		}
	}

return cmed;
}


int flip_spin(int *red,int dim,int *spin,float *delta_E,float J,float B,float *p_l,float *p_i){ 
	int i,ind_l,ind_i,cmed;
	float p;
	i=(int)(round(aleatorio()*((float)(dim*dim))));
	*spin=*(red+i);
		
	cmed=campo_med(red,dim,i);
	
	*delta_E=2*(*spin)*(B+J*cmed); 	
	
	ind_l=(*spin+1)/2;								//indice de la probabilidad libre
	ind_i=((*spin)*cmed+4)/2;							//indeice de la probabilidad de interaccion
	if(*delta_E<0.0){	
		*(red+i)=(*spin)*(-1);}
	else{	
		p=(*(p_i+ind_i))*(*(p_l+ind_l));
		if(aleatorio()<p){
			*(red+i)=(*spin)*(-1);
		}
	}

return 0;	
}	

int flip_libre(int *red,int dim,int *spin,float *delta_E,float B,float *p_l,int *indice){ 
	int i,ind_l,ind_i,cmed;
	float p;
	i=(int)(round(aleatorio()*((float)(dim*dim))));
	*spin=*(red+i);
		
	*indice=i;
	*delta_E=2*(*spin)*B; 	
	
	ind_l=(*spin+1)/2;								//indice de la probabilidad libre
	
	
		p=*(p_l+ind_l);
//		printf("%lf %lf\n",p, *delta_E, ind_l);
		if(aleatorio()<p){
			*(red+i)=(*spin)*(-1);
		}
		else{
			*delta_E=0.0;
		}
return 0;	
}	

int poblar(int *red, float p, int dim){
	int i;
	for(i=0;i<dim*dim;i++){
		*(red+i)=-1;
		if(aleatorio()<p){
			*(red+i)=1;		
		}	
	}
return 0;
}


float M(int *red, int dim){
	int i;
	float m;

	m=0;	

	for(i=0;i<dim*dim;i++){
		m=m+(*(red+i));
	}
return m/((float)(dim*dim));
}

float energia_sin_int(int *red, int dim, float B){
	int i;
	float u;

	u=0;	

	for(i=0;i<dim*dim;i++){
		u=u-*(red+i);
	}
return B*u;
}

int energia_con_int(int *red, int dim, float B,float J){
	int i;
	float u;

	u=0;	

	for(i=0;i<dim*dim;i++){
		u=u-*(red+i);
	}
return B*u;
}




int problema2a(int dim,int N){
	int i, j, k, *red,*spin, N_prom,*indice;
	float p,*u, aceptacion,*m,*m_cuad,*H,*H_cuad,*delta_E,*p_l,B0,B,delta_T,T,a,b;
	
	FILE *fp= fopen("1a_dim10", "w");
	FILE *fc= fopen("correlacion2a", "w");
	FILE *findi= fopen("indices", "a");	
	

	delta_T=0.1;	
	aceptacion=0;
	N_prom=dim*dim;
		
	red=(int*)malloc(dim*dim*sizeof(int));
	delta_E=(float*)malloc(sizeof(float));
	spin=(int*)malloc(sizeof(int));
	indice=(int*)malloc(sizeof(int));
	
	m=(float*)malloc(sizeof(float));
	H=(float*)malloc(sizeof(float));
	m_cuad=(float*)malloc(sizeof(float));
	H_cuad=(float*)malloc(sizeof(float));
	p_l=(float*)malloc(2*sizeof(float));

	poblar(red, 0.5, dim);

	
	B=0.5;
	p_libre(B,p_l);
	

	
	for(B=0.0;B<10.0;B=B+0.1){
		fprintf(fc, "%f",B);	
		//B=B0/T;
		p_libre(B,p_l);
		for(i=0;i<100000;i++){
			flip_libre(red,dim,spin,delta_E,B,p_l,indice);    //termzlizacion
			fprintf(fc, " %f",M(red,dim));
			}
			fprintf(fc, "\n");

		*m=0.0;
		*H=0.0;
		for(i=0;i<N_prom;i++){
			flip_libre(red,dim,spin,delta_E,B,p_l,indice);
			a=energia_sin_int(red,dim, B);
			b=M(red,dim);
			*H=*H+a;
			*H_cuad=*H_cuad+a*a;
			*m=*m+b;
			*m_cuad=*m_cuad+b*b;}
		*H=*H/((float)N_prom);
		*H_cuad=*H_cuad/((float)N_prom);
		*m=*m/((float)N_prom);
		*m_cuad=*m_cuad/((float)N_prom);
		//imprimir(red,dim);
		fprintf(fp, "%f %f %f %f %f %f\n", 1.0/B,B, *H, *m,*H_cuad,*m_cuad);
		
	}

free(red);
free(delta_E);
free(spin);
free(m);
free(H);
free(p_l);
fclose(fp);
fclose(fc);
fclose(findi);
return 0;
}



int problema2b(int N){
	int i, j, k, *red, dim,paso,*spin;
	float p,*probas,J,*u, acum,*m,*H,*delta_E,*p_i,*p_l,B;
		
	dim=10;
	paso=0;	
	acum=0;
		
	red=(int*)malloc(dim*dim*sizeof(int));
	probas=(float*)malloc(5*sizeof(float));
	delta_E=(float*)malloc(sizeof(float));
	spin=(int*)malloc(sizeof(int));

	
	m=(float*)malloc(N*sizeof(float));
	H=(float*)malloc(N*sizeof(float));
	p_i=(float*)malloc(4*sizeof(float));
	p_l=(float*)malloc(2*sizeof(float));

	poblar(red, 0.5, dim);

	
	J=0.0;
	B=1.0;
	
	p_libre(B,p_i);
	p_int(J,p_l);
	
	
	*H=energia_sin_int(red,dim, B);
	*m=M(red,dim);
	for(i=0;i<N-1;i++){
		flip_spin(red,dim,spin,delta_E,J,B,p_l,p_i);
		*(H+i+1)=*(H+i)+*delta_E;
		*(m+i+1)=*(m+i)+2*(*spin);
	}
	
	


	FILE *fp= fopen("a", "w");


free(red);
fclose(fp);
return 0;
}









int imprimir(int *red, int dim){
	int i,j;
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			if(*(red+dim*i+j)<10){printf("%d   ",*(red+dim*i+j));}
			else{printf("%d  ",*(red+dim*i+j));}
		}	
	printf("\n");
	}
printf("\n");
return 0;
}

