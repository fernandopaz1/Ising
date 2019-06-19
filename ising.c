#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>



int rand(void);
double aleatorio(void);
int poblar(int *red, double p, int dim);
int imprimir(int *red, int dim);
int masa(int *red, int dim,int *etiqueta_percolante, int *mass);
int problema2a(int dim, int N);
int problema2b(int dim,int N);
int problema2d(int dim,int N);
int problema2e(int dim,int N);
int campo_med(int *red,int dim,int i,int j);
int campo_med_2do_vecinos(int *red,int dim,int i,int j);
double M(int *red, int dim);
double energia_con_int(int *red, int dim, double B,double J);
double energia_sin_int(int *red, int dim, double B);
int p_libre(double B,double *p_l);
int p_int(double J,double *p_i);
int flip_spin(int *red,int dim,int *spin,double *delta_E,double J,double B,double *p_l,double *p_i,double *m,double *H);
int flip_antiferro(int *red,int dim,int *spin,double *delta_E,double J,double B,double *p_l,double *p_i,double *m,double *H);
int flip_frustracion(int *red,int dim,int *spin,double *delta_E,double J,double mu,double *p_l,double *p_i,double *m,double *H);
int flip_libre(int *red,int dim,int *spin,double *delta_E,double B,double *p_l,int *indice);
double mean(double *a,int inicio,int lenght);
double std_dev(double *a, int inicio,int lenght, double mean_sqrt);
double correl(double *a,int inicio, int paso, int lenght, int porcion,double media, double desvio);
int correlacion(int dim,int iteraciones);


int main(int argc,char *argv[]){
	int iteraciones,dim;


	double total_time;
	clock_t start, end;
	start = clock();

	sscanf(argv[1], "%d", & dim);
	sscanf(argv[2], "%d", & iteraciones);
	
	srand(time(NULL)); 



//	problema2a(dim,iteraciones);
//	problema2b(dim,iteraciones);
//	problema2d(dim,iteraciones);
	problema2e(dim,iteraciones);
//	correlacion(dim,iteraciones);

	end = clock();								//time count stops 


	total_time = ((double) (((double)(end - start)) / (double)CLOCKS_PER_SEC)/60.0);	
	printf("\nEl tiempo (minutos) requerido es:  %lf \n", total_time);		//calulate total time

	
	return 0;
}





double aleatorio(void){
	double a;
	a=((double)rand())/((double)RAND_MAX);
return a;
}



int p_int(double J,double *p_i){                           //Cuando B=0   delta_E es cmed*J y cmed va desde 8.0 a -8.0
	
	*(p_i)=exp(8.0*J);
	*(p_i+1)=exp(4.0*J);
	*(p_i+2)=1.0;
	*(p_i+3)=exp(-4.0*J);
	*(p_i+4)=exp(-8.0*J);
return 0;
}



int p_libre(double B,double *p_l){
	int i;
	for(i=0;i<2;i++){
		*(p_l+i)=exp(2.0*(1.0-2.0*((double)i))*B);
	}
return 0;
}

int campo_med(int *red,int dim,int i,int j){
	int cmed,up,down,left,right;
	
	up=(i-1+dim)%dim;
	down=(i+1+dim)%dim;
	left=(j-1+dim)%dim;
	right=(j+1+dim)%dim;
	
	cmed=*(red+dim*up+j)+*(red+dim*down+j)+*(red+dim*i+left)+*(red+dim*i+right);
	

return cmed;
}

int campo_med_2do_vecinos(int *red,int dim,int i,int j){
	int cmed,up,down,left,right;
	
	up=(i-1+dim)%dim;
	down=(i+1+dim)%dim;
	left=(j-1+dim)%dim;
	right=(j+1+dim)%dim;
	
	cmed=*(red+dim*up+left)+*(red+dim*up+right)+*(red+dim*down+left)+*(red+dim*down+right);

return cmed;
}



int flip_spin(int *red,int dim,int *spin,double *delta_E,double J,double B,double *p_l,double *p_i,double *m,double *H){ 
	int i,j,ind_l,ind_i,cmed;
	double p;
	i=(int)(aleatorio()*dim);
	j=(int)(aleatorio()*dim);
	*spin=*(red+dim*i+j);
		
	cmed=campo_med(red,dim,i,j);
	
	*delta_E=2.0*((double)(*spin))*(B+J*((double)cmed)); 	
	
		

	ind_l=(*spin+1)/2;								//indice de la probabilidad libre
	ind_i=((*spin)*cmed+4)/2;							//indice de la probabilidad de interaccion
	p=(*(p_i+ind_i))*(*(p_l+ind_l));						//eval la porbabilidad
	if(aleatorio()<=p){
		*(red+dim*i+j)=(*spin)*(-1);
		*H+=*delta_E;                                                         //calcula el H en cada flipeo
		*m-=((*spin)*2);							//Calcula M*dim*dim, guarda luego hay que normalizar
	}

return 0;	
}	

int flip_libre(int *red,int dim,int *spin,double *delta_E,double B,double *p_l,int *indice){ 
	int i,ind_l,ind_i,cmed;
	double p;
	i=(int)(round(aleatorio()*((double)(dim*dim))));
	*spin=*(red+i);
		
	*indice=i;
	*delta_E=2*(*spin)*B; 	
	
	ind_l=(*spin+1)/2;								//indice de la probabilidad libre
	
	
	p=*(p_l+ind_l);
	
//	printf("%lf %lf\n",p, *delta_E, ind_l);
	if(aleatorio()<p){
		*(red+i)=(*spin)*(-1);

	}
return 0;	
}	

int flip_antiferro(int *red,int dim,int *spin,double *delta_E,double J,double mu,double *p_l,double *p_i,double *m,double *H){ 
	int i,j,ind_l,ind_i,cmed;
	double p;
	i=(int)(aleatorio()*dim);
	j=(int)(aleatorio()*dim);
	*spin=*(red+dim*i+j);
		
	cmed=campo_med(red,dim,i,j);
	
	*delta_E=2.0*((double)(*spin))*J*(mu+((double)cmed)); 	
	
		

	ind_l=(*spin+1)/2;								//indice de la probabilidad libre
	ind_i=((*spin)*cmed+4)/2;							//indice de la probabilidad de interaccion
	p=(*(p_i+ind_i))*(*(p_l+ind_l));						//eval la porbabilidad
	if(aleatorio()<=p){
		*(red+dim*i+j)=(*spin)*(-1);
		*H+=*delta_E;                                                         //calcula el H en cada flipeo
		*m-=((*spin)*2);							//Calcula M*dim*dim, guarda luego hay que normalizar
	}

return 0;	
}	


int flip_frustracion(int *red,int dim,int *spin,double *delta_E,double J,double mu,double *p_l,double *p_i,double *m,double *H){ 
	int i,j,ind_l,ind_i,cmed,cmed2do,ind_i_antiferro;
	double p;
	i=(int)(aleatorio()*dim);
	j=(int)(aleatorio()*dim);
	*spin=*(red+dim*i+j);
		
	cmed=campo_med(red,dim,i,j);
	cmed2do=campo_med_2do_vecinos(red,dim,i,j);
	*delta_E=2.0*((double)(*spin))*J*(mu+((double)(cmed-cmed2do))); 	
	
		

	ind_l=(*spin+1)/2;								//indice de la probabilidad libre
	ind_i=((*spin)*cmed+4)/2;							//indice de la probabilidad de interaccion
	ind_i_antiferro=(-(*spin)*cmed2do+4)/2;						//indice de la probabilidad a segundo vecino
	p=(*(p_i+ind_i))*(*(p_l+ind_l))*(*(p_i+ind_i_antiferro));			//eval la porbabilidad
	if(aleatorio()<=p){
		*(red+dim*i+j)=(*spin)*(-1);
		*H+=*delta_E;                                                         //calcula el H en cada flipeo
		*m-=((*spin)*2);							//Calcula M*dim*dim, guarda luego hay que normalizar
	}

return 0;	
}	

int poblar(int *red, double p, int dim){
	int i;
	for(i=0;i<dim*dim;i++){
		if(aleatorio()<p){
			*(red+i)=1;		
		}
		else{*(red+i)=-1;}	
	}
return 0;
}


double M(int *red, int dim){
	int i;
	double m;

	m=0.0;	

	for(i=0;i<dim*dim;i++){
		m+=(double)(*(red+i));
	}
return ((double)m);
}

double energia_sin_int(int *red, int dim, double B){
	int i;
	double u;

	u=0.0;	

	for(i=0;i<dim*dim;i++){
		u-=*(red+i);
	}
return B*u;
}

double energia_con_int(int *red, int dim, double B,double J){
	int i,j;
	double u;

	u=0.0;	

	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			u-=((double)(*(red+dim*i+j))*((double)(B+J*campo_med(red,dim,i, j))));
		}
	}
	//printf("%lf\n",u);
return u;
}


double mean_cuad(double *a,int inicio,int lenght){
	int i;  
	double sum=0;  
	for(i=0;i<lenght;i++){ 
		sum+=(double)(*(a+inicio+i));}
	sum=sum*sum/((double)(lenght*lenght));
	printf("\n media cuadrada %lf",sum);
return sum;}

double std_dev(double *a, int inicio,int lenght, double mean_sqrt){
	int i;  
	double sum=0;  
	for(i=0;i<lenght;i++){ 
		sum+=(double)(*(a+inicio+i))*(*(a+inicio+i));}
	sum=(sum/((double)lenght));
	printf("\n Desvio %lf",(sum/((double)lenght)));
	printf("\n Desvio2 %lf",sum-mean_sqrt);
	printf("\n");

return (sum-mean_sqrt);}

double correl(double *a,int inicio, int paso, int lenght, int porcion,double media_cuad, double desvio){
	int j;  
	double sum=0.0;
	for(j=0;j<(porcion-paso);j++){
		sum=sum+(double)(((*(a+inicio+j))*(*(a+inicio+j+paso))));
	}
	sum=sum/((double)(porcion-paso)); 
	//	printf("\n Desvio %lf",sum);
return (sum-media_cuad)/(desvio);
}


int correlacion(int dim,int N){
	int i, j, k, *red,*spin, N_prom,*indice,N_div,N_total;
	double p,*u,*m_cuad,*H_cuad,*delta_E,*p_l,*p_i,B0,B,J,delta_T,T,a,b,*m_acum,*H_acum,porcion,*rho_m,*rho_h;
	double 	*m,*H,media_cuad,desvio;
	
	FILE *fterm= fopen("termalizacion2b", "w");
	

	N_prom=dim*dim;
		
	red=(int*)malloc(dim*dim*sizeof(int));
	delta_E=(double*)malloc(sizeof(double));
	spin=(int*)malloc(sizeof(int));
	indice=(int*)malloc(sizeof(int));
	N_div=1000000;	
	N_total=10*N_div;
	porcion=N_div/20;

	m=(double*)malloc(sizeof(double));
	H=(double*)malloc(sizeof(double));
	rho_m=(double*)malloc(porcion*sizeof(double));
	rho_h=(double*)malloc(porcion*sizeof(double));
	m_acum=(double*)malloc(N_total*sizeof(double));
	H_acum=(double*)malloc(N_total*sizeof(double));
	m_cuad=(double*)malloc(sizeof(double));
	H_cuad=(double*)malloc(sizeof(double));
	p_l=(double*)malloc(2*sizeof(double));
	p_i=(double*)malloc(5*sizeof(double));

	poblar(red, 0.5, dim);


	
	
	
	J=0.1;
	B=0.0;
	
	p_libre(0,p_l);
	p_int(J,p_i);
	
	J=0.1;
	a=((double)dim*dim);
	double a_cuad=a*a;
	*H=energia_con_int(red,dim, B,J);
	*m=M(red,dim);

	for(i=0;i<dim*dim;i++){
			flip_spin(red,dim,spin,delta_E,J, 0.0, p_l,p_i,m,H);}	


	for(J=0.1;J<0.6;J=J+0.1){

		fprintf(fterm, "%lf",J);	
		p_int(J,p_i);
		
		for(i=0;i<10*N_div;i++){
			flip_spin(red,dim,spin,delta_E,J, 0.0, p_l,p_i,m,H);    //termzlizacion
			//printf("%lf\n",*(m_acum+i)=M(red,dim));
			*(m_acum+i)=*m/a;
			*(H_acum+i)=*H;
		//	fprintf(fterm, " %lf %lf",*m,*H);
		//	fprintf(fterm, " %lf %lf",M(red,dim),energia_con_int(red,dim, B,J));
			}
		for(i=0;i<10;i++){
			media_cuad=mean_cuad(m_acum, i*N_div,N_div)/a_cuad;
			desvio=std_dev(m_acum,i*N_div, N_div,media_cuad)/a_cuad;
			for(j=0;j<porcion;j++){
				*(rho_m+j)=*(rho_m+j)+correl(m_acum,i*N_div, j, N_div,porcion, media_cuad, desvio)/10.0;
			}			
		}		
		i=0;
		//media_cuad=mean_cuad(m_acum, i*N_div,N_div);
		//desvio=std_dev(m_acum,i*N_div, N_div,media_cuad);
		printf("%lf",1.0);
		for(j=0;j<porcion;j++){
				fprintf(fterm, " %lf",*(rho_m+j));
	//			fprintf(fterm, " %lf",correl(m_acum,i*N_div, j, N_div,porcion, media_cuad, desvio));
			}
			
		fprintf(fterm,"\n");

	}

fclose(fterm);
return 0;
}


int problema2a(int dim,int N){
	int i, j, k, *red,*spin, N_prom,*indice;
	double p,*u, aceptacion,*m,*m_cuad,*H,*H_cuad,*delta_E,*p_l,B0,B,delta_T,T,a,b;
	
	FILE *fp= fopen("1a_dim10", "w");
	FILE *findi= fopen("indices", "a");	
	

	delta_T=0.1;	
	aceptacion=0;
	N_prom=dim*dim;
		
	red=(int*)malloc(dim*dim*sizeof(int));
	delta_E=(double*)malloc(sizeof(double));
	spin=(int*)malloc(sizeof(int));
	indice=(int*)malloc(sizeof(int));
	
	m=(double*)malloc(sizeof(double));
	H=(double*)malloc(sizeof(double));
	m_cuad=(double*)malloc(sizeof(double));
	H_cuad=(double*)malloc(sizeof(double));
	p_l=(double*)malloc(2*sizeof(double));

	poblar(red, 0.5, dim);

	
	B=0.5;
	p_libre(B,p_l);
	

	
	for(B=0.0;B<10.0;B=B+0.1){
		//fprintf(fp, "%lf",B);	
		//B=B0/T;
		p_libre(B,p_l);
		for(i=0;i<100000;i++){
			flip_libre(red,dim,spin,delta_E,B,p_l,indice);    //termzlizacion
			fprintf(fp, " %lf",M(red,dim));
			}
			fprintf(fp, "\n");

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
		*H=*H/((double)N_prom);
		*H_cuad=*H_cuad/((double)N_prom);
		*m=*m/((double)N_prom);
		*m_cuad=*m_cuad/((double)N_prom);
		//imprimir(red,dim);
		fprintf(fp, "%lf %lf %lf %lf %lf %lf\n", 1.0/B,B, *H, *m,*H_cuad,*m_cuad);
		
	}

free(red);
free(delta_E);
free(spin);
free(m);
free(H);
free(p_l);
fclose(fp);
fclose(findi);
return 0;
}



int problema2b(int dim,int N){
	int i, j, k, *red,paso,*spin,N_prom;
	double p,J,*u, acum,*delta_E,*p_i,*p_l,B;
	double	*m,*H,*H_cuad,*m_cuad,a,b,*m_acum,*H_acum,*H_cuad_acum,*m_cuad_acum;
	
	
	paso=0;	
	acum=0;
	N_prom=10000;
		

	FILE *fp= fopen("2b", "w");

	red=(int*)malloc(dim*dim*sizeof(int));
	delta_E=(double*)malloc(sizeof(double));
	spin=(int*)malloc(sizeof(int));

	
	m=(double*)malloc(sizeof(double));
	H=(double*)malloc(sizeof(double));
	m_cuad=(double*)malloc(sizeof(double));
	H_cuad=(double*)malloc(sizeof(double));

	m_acum=(double*)malloc(sizeof(double));
	H_acum=(double*)malloc(sizeof(double));
	m_cuad_acum=(double*)malloc(sizeof(double));
	H_cuad_acum=(double*)malloc(sizeof(double));

	p_i=(double*)malloc(5*sizeof(double));
	p_l=(double*)malloc(2*sizeof(double));

	poblar(red, 0.5, dim);         //puebla la red

	
	J=0.2;	
	double B0=0.0;         
	B=B0/J;
	p_libre(B0,p_l);
	p_int(J,p_i);
	
	
	*H=energia_con_int(red,dim, B,J);
	*m=M(red,dim);
	while(J<1.0){
		B=-B0/J;
		p_libre(B,p_l);		//actualiza las probabilidades
		p_int(J,p_i);
		

		*m_acum=0.0;		//seteo contadores a cero
		*H_acum=0.0;
		*m_cuad_acum=0.0;
		*H_cuad_acum=0.0;
		a=0.0;
		b=0.0;
		
		for(i=0;i<N_prom;i++){

			for(k=0;k<5*dim*dim;k++){
				flip_spin(red,dim,spin,delta_E,J, B, p_l,p_i,m,H);    //termzlizacion
			}

			a=*H;
			b=*m;			// en cada flipeo de spin calcula M*dim*dim y H usando *spin y *delta_E

			*H_acum=*H_acum+a;
			*H_cuad_acum=*H_cuad_acum+a*a;
			*m_acum=*m_acum+b;
			*m_cuad_acum=*m_cuad_acum+b*b;
		}

		*H_acum=*H_acum/((double)N_prom);
		*H_cuad_acum=*H_cuad_acum/((double)N_prom);
		*m_acum=*m_acum/((double)(N_prom*dim*dim));			//Divide M por N_prom y por dim*dim
		*m_cuad_acum=*m_cuad_acum/((double)(N_prom*dim*dim));		//Divide M**2 por N_prom y por dim*dim*dim*dim
		*m_cuad_acum=*m_cuad_acum/((double)(dim*dim));			//lo hago en dos tandas porque dim*dim*dim*dim es muy grande
	

		fprintf(fp, "%lf %lf %lf %lf %lf %lf\n", 1.0/J,J, *H_acum, *m_acum,*H_cuad_acum,*m_cuad_acum);

		printf("%lf\n",J);
		if(J<0.4 || J>0.8){			//paso en J comun
			J=J+0.005;
		}
		else{J=J+0.0005;}			//dismuniyo en la zona critica
	}
	


	

free(red);
free(delta_E);
free(spin);
free(m);
free(H);
free(m_cuad);
free(H_cuad);
free(p_l);
free(p_i);
fclose(fp); 
return 0;
}






int problema2d(int dim,int N){
	int i, j, k, *red,paso,*spin,N_prom;
	double p,J,T,*u, acum,*delta_E,*p_i,*p_l,B;
	double	*m,*H,*H_cuad,*m_cuad,a,b,*m_acum,*H_acum,*H_cuad_acum,*m_cuad_acum,*B_tira;
	
	
	paso=0;	
	acum=0;
	N_prom=1000;
		

	FILE *fp= fopen("2d", "w");

	red=(int*)malloc(dim*dim*sizeof(int));
	delta_E=(double*)malloc(sizeof(double));
	spin=(int*)malloc(sizeof(int));

	
	m=(double*)malloc(sizeof(double));
	H=(double*)malloc(sizeof(double));
	m_cuad=(double*)malloc(sizeof(double));
	H_cuad=(double*)malloc(sizeof(double));

	m_acum=(double*)malloc(sizeof(double));
	H_acum=(double*)malloc(sizeof(double));
	m_cuad_acum=(double*)malloc(sizeof(double));
	H_cuad_acum=(double*)malloc(sizeof(double));

	p_i=(double*)malloc(5*sizeof(double));
	p_l=(double*)malloc(2*sizeof(double));

	
	
	int pasos=100;
	double B_max=0.50;
	double incremento=B_max/((double)pasos);
	B_tira=(double*)malloc(pasos*5*sizeof(double));

	poblar(red, 0.5, dim);         //puebla la red


	//double J0=-0.5;	
	double B0=1.0; 
	//T=1.0;        
	
	J=0.2;
	B=B0/J;
	p_libre(B0,p_l);
	p_int(J,p_i);
	

	double mu=0.0;	

	incremento=0.005;
	*H=energia_con_int(red,dim, B,J);
	*m=M(red,dim);
	for(J=-0.2;J>-1.0;J-=incremento){
		B=fabs(B0/J);
		p_libre(B,p_l);		//actualiza las probabilidades
		p_int(J,p_i);
		

		*m_acum=0.0;		//seteo contadores a cero
		*H_acum=0.0;
		*m_cuad_acum=0.0;
		*H_cuad_acum=0.0;
		a=0.0;
		b=0.0;
		
		for(i=0;i<N_prom;i++){

			for(k=0;k<5*dim*dim;k++){
				flip_antiferro(red,dim,spin,delta_E,J, mu, p_l,p_i,m,H);    //termzlizacion
			}

			a=*H;
			b=*m;			// en cada flipeo de spin calcula M*dim*dim y H usando *spin y *delta_E

			*H_acum=*H_acum+a;
			*H_cuad_acum=*H_cuad_acum+a*a;
			*m_acum=*m_acum+b;
			*m_cuad_acum=*m_cuad_acum+b*b;
		}

		*H_acum=*H_acum/((double)N_prom);
		*H_cuad_acum=*H_cuad_acum/((double)N_prom);
		*m_acum=*m_acum/((double)(N_prom*dim*dim));			//Divide M por N_prom y por dim*dim
		*m_cuad_acum=*m_cuad_acum/((double)(N_prom*dim*dim));		//Divide M**2 por N_prom y por dim*dim*dim*dim
		*m_cuad_acum=*m_cuad_acum/((double)(dim*dim));			//lo hago en dos tandas porque dim*dim*dim*dim es muy grande
	

		fprintf(fp, "%lf %lf %lf %lf %lf %lf\n", -1.0/J,J, *H_acum, *m_acum,*H_cuad_acum,*m_cuad_acum);

		//printf("%lf\n",J);
		printf("%lf\n",J);

		if(J<-0.4 || J>-0.8){		//dismuniyo en la zona critica
			J=J-0.1*incremento;
		}
		
	}
	



free(red);
free(delta_E);
free(spin);
free(m);
free(H);
free(B_tira);
free(m_cuad);
free(H_cuad);
free(p_l);
free(p_i);
fclose(fp); 
return 0;
}




int problema2e(int dim,int N){
	int i, j, k, *red,paso,*spin,N_prom;
	double p,J,T,*u, acum,*delta_E,*p_i,*p_l,B;
	double	*m,*H,*H_cuad,*m_cuad,a,b,*m_acum,*H_acum,*H_cuad_acum,*m_cuad_acum,*B_tira;
	
	
	paso=0;	
	acum=0;
	N_prom=1000;
		

	FILE *fp= fopen("2e", "w");

	red=(int*)malloc(dim*dim*sizeof(int));
	delta_E=(double*)malloc(sizeof(double));
	spin=(int*)malloc(sizeof(int));

	
	m=(double*)malloc(sizeof(double));
	H=(double*)malloc(sizeof(double));
	m_cuad=(double*)malloc(sizeof(double));
	H_cuad=(double*)malloc(sizeof(double));

	m_acum=(double*)malloc(sizeof(double));
	H_acum=(double*)malloc(sizeof(double));
	m_cuad_acum=(double*)malloc(sizeof(double));
	H_cuad_acum=(double*)malloc(sizeof(double));

	p_i=(double*)malloc(5*sizeof(double));
	p_l=(double*)malloc(2*sizeof(double));

	
	
	int pasos=20;

	poblar(red, 0.5, dim);         //puebla la red



	double B0=1.0;

 
	double mu=0.0; 


	J=0.2;
	B=mu*J;
	p_libre(B,p_l);
	p_int(J,p_i);
	

	pasos=0.01;
	*H=energia_con_int(red,dim, B,J);
	*m=M(red,dim);
	for(J=0.2;J<3.0;J+=0.01){
		B=J*mu;
		p_libre(B,p_l);		//actualiza las probabilidades
		p_int(J,p_i);
		

		*m_acum=0.0;		//seteo contadores a cero
		*H_acum=0.0;
		*m_cuad_acum=0.0;
		*H_cuad_acum=0.0;
		a=0.0;
		b=0.0;
		
		for(i=0;i<N_prom;i++){

			for(k=0;k<5*dim*dim;k++){
				flip_frustracion(red,dim,spin,delta_E,J, mu, p_l,p_i,m,H);    //termzlizacion
			}

			a=*H;
			b=*m;			// en cada flipeo de spin calcula M*dim*dim y H usando *spin y *delta_E

			*H_acum=*H_acum+a;
			*H_cuad_acum=*H_cuad_acum+a*a;
			*m_acum=*m_acum+b;
			*m_cuad_acum=*m_cuad_acum+b*b;
		}

		*H_acum=*H_acum/((double)N_prom);
		*H_cuad_acum=*H_cuad_acum/((double)N_prom);
		*m_acum=*m_acum/((double)(N_prom*dim*dim));			//Divide M por N_prom y por dim*dim
		*m_cuad_acum=*m_cuad_acum/((double)(N_prom*dim*dim));		//Divide M**2 por N_prom y por dim*dim*dim*dim
		*m_cuad_acum=*m_cuad_acum/((double)(dim*dim));			//lo hago en dos tandas porque dim*dim*dim*dim es muy grande
	

		fprintf(fp, "%lf %lf %lf %lf %lf %lf\n", 1.0/J,J, *H_acum, *m_acum,*H_cuad_acum,*m_cuad_acum);

		printf("%lf\n",J);
		

	/*	if(J>0.4 || J<0.8){		//dismuniyo en la zona critica
			J=J-0.1*incremento;
		}
	*/	
	}
	



free(red);
free(delta_E);
free(spin);
free(m);
free(H);
free(m_cuad);
free(H_cuad);
free(p_l);
free(p_i);
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

