#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>



int rand(void);
float aleatorio(void);
float powf(float x, float y);
int poblar(int *red, float p, int dim);
//int imprimir(int *red, int dim);
int masa(int *red, int dim,int *etiqueta_percolante, int *mass);
int problema2(int N);
int campo_med(int *red,int dim,int i);
int p_libre(float J,float *p_l);
int p_int(float B,float *p_i);
int flip_spin(int *red,int dim,int *spin,float *delta_E,float J,float B,float *p_l,float *p_i);

int main(int argc,char *argv[]){
	int iteraciones;


	double total_time;
	clock_t start, end;
	start = clock();

	
	sscanf(argv[1], "%d", & iteraciones);
	srand(time(NULL)); 



	problema2(iteraciones);

	end = clock();								//time count stops 


	total_time = ((double) (((end - start)) / CLOCKS_PER_SEC)/60)/60;	
	printf("\nEl tiempo (horas) requerido es:  %f \n", total_time);		//calulate total time

	
	return 0;
}





float aleatorio(void){
	float a;
	a=((float)rand())/((float)RAND_MAX);
return a;
}



int p_libre(float J,float *p_l){
	int i;
	for(i=0;i<5;i++){
		*p_l=exp(-8*(1.0-0.5*((float)(i)))*J);
	}
return 0;
}

int p_int(float B,float *p_i){
	int i;
	for(i=0;i<2;i++){
		*p_i=exp(-2*(1.0-0.5*((float)(i)))*B);
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
	i=round(aleatorio()*((float)(dim*dim)));
	*spin=*(red+i);
	cmed=campo_med(red,dim,i);
	
	*delta_E=2*(*spin)*(B+J*cmed); 	
	
	ind_l=(*spin+1)/2;								//indice de la probabilidad libre
	ind_i=((*spin)*cmed+4)/2;							//indeice de la probabilidad de interaccion
	if(*delta_E<0){	
		*(red+i)=(-1)**(red+i);}
	else{	
		p=(*(p_i+ind_i))*(*(p_l+ind_l));
		if(aleatorio()<p){
			*(red+i)=(-1)**(red+i);
		}
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


int M(int *red, int dim){
	int i;
	float m,u;

	m=0;	

	for(i=0;i<dim*dim;i++){
		m=m+(*(red+i));
	}
return m/(dim*dim);
}

int energia_sin_int(int *red, int dim, float J){
	int i;
	float u;

	u=0;	

	for(i=0;i<dim*dim;i++){
		u=u+J*(*(red+i));
	}
return u;
}




int problema2(int N){
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
	p_i=(float*)malloc(5*sizeof(float));
	p_l=(float*)malloc(2*sizeof(float));

	poblar(red, 0.5, dim);

	
	J=2.0;
	B=0.0;
	
	p_libre(B,p_i);
	p_int(J,p_l);
	
	
	*H=energia_sin_int(red,dim, J);
	*m=M(red,dim);
	for(i=0;i<N-1;i++){
		flip_spin(red,dim,spin,delta_E,J,B,p_l,p_i);
		*(H+i+1)=*(H+i)+*delta_E;
	}
	
	


	FILE *fp= fopen("pc_nu", "w");


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
return 0;
}

