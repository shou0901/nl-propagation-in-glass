#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_fft_complex.h>
#define N 8192 // FFT
#define M 100  // dz = L/M
#define PI 3.1415926535897932384626433832795

// gcc pulseKonishi.c -L/usr/local/lib -lgsl -lgslcblas

double getEnergy(double *Cdata, double dt);

double n2 = 3.5e-20; // [m^2/W] Fused Silica Boyd p212
double c = 2.99792458e8; // [m/s]
double wavelength = 470.0e-9; // [m]
double epsilon0 = 8.854e-12; // permittivity of free space
double n0 = 1.7771;
double tau = 45.0e-15; // FWHM of pulse duration[s]
double beam_w = 50e-6; // beam radius [m]
double E = 1.0; // pulse energy [J]
double L = 500.0e-6; // [m]thickness
double k2 = 122.75e-27; // [s^2/m] refractiveindex.info Fused Silica Malitson 1965
double k3 = 0.0*27.47e-42; // [s^2/m] newport

int main()
{
	// FFT http://jajagacchi.web.fc2.com/dft.html
int i,j,k;
double temp,phase;
FILE *fp;
char str[100];

double omega0 = 2.0*PI*c/wavelength; // frequency [rad/s]
//double chi3 = 3.78e-24; //m^2/V^2
double gamma = 2.0*n0*epsilon0*n2*omega0;
double I0; // peak intensity
double F0;// peak electric field amplitude

double t,omega;
double start_t = -100.0*tau;
double end_t = 100.0*tau;
double dt = (end_t-start_t)/(double)N;
double dz = L/(double)M;
double domega = 2.0*PI/(double)N/dt;
double start_omega = omega0-(double)(N/2)*domega;

double *Cdata = (double*)malloc(sizeof(double)*2*N);
double *temp1 = (double*)malloc(sizeof(double)*2*N);
double *temp2 = (double*)malloc(sizeof(double)*2*N);
gsl_fft_complex_wavetable *wavetable = gsl_fft_complex_wavetable_alloc(N);
gsl_fft_complex_workspace *workspace = gsl_fft_complex_workspace_alloc(N);

for(k=0;k<1;k++)
{
	E = (double)(25-k)*1.0e-6;
	I0 = 2.0*E/PI/beam_w/beam_w*sqrt(4.0*log(2.0)/PI)/tau; // peak intensity
	F0 = sqrt(2.0*I0/epsilon0/c); // peak electric field amplitude

for(i=0;i<N;i++)
{
	t = start_t+dt*(double)i;
	Cdata[2*i] = F0*exp(-2.0*log(2.0)*t*t/tau/tau);
	Cdata[2*i+1] = 0.0;
	phase = (start_omega-omega0)*t;
	temp1[2*i] = (Cdata[2*i]*cos(phase)+Cdata[2*i+1]*sin(phase))*dt/sqrt(2.0*PI);
	temp1[2*i+1] = (-Cdata[2*i]*sin(phase)+Cdata[2*i+1]*cos(phase))*dt/sqrt(2.0*PI);
}

sprintf(str,"beforeKonishi%d.dat",k);
fp = fopen(str,"w");
for(i=0;i<N;i++)
{
	t = start_t+dt*(double)i;
	fprintf(fp,"%e %e\n",t*1.0e15,c*epsilon0*(Cdata[2*i]*Cdata[2*i]+Cdata[2*i+1]*Cdata[2*i+1])/2.0*1.0e-4);
}
fclose(fp);

gsl_fft_complex_forward(temp1,1,N,wavetable,workspace);
sprintf(str,"beforeKonishispectrum%d.dat",k);
fp = fopen(str,"w");
for(j=0;j<N;j++)
{
	omega = start_omega + domega*(double)j;
	fprintf(fp,"%e %e\n",omega*1.0e-15,temp1[2*j]*temp1[2*j]+temp1[2*j+1]*temp1[2*j+1]);
}
fclose(fp);

//Split-Step Fourier Method

for(i=0;i<M;i++)
{
	// temp copy
	for(j=0;j<N;j++)
	{
		temp1[2*j] = Cdata[2*j];
		temp1[2*j+1] = Cdata[2*j+1];
		temp2[2*j] = Cdata[2*j];
		temp2[2*j+1] = Cdata[2*j+1];
	}

	//delta_z/2 proceed
	for(j=0;j<N;j++)
	{
		temp = temp2[2*j]*temp2[2*j]+temp2[2*j+1]*temp2[2*j+1];
		t = start_t+dt*(double)j;
		phase = 0.5*gamma*temp*dz+(start_omega-omega0)*t;
		Cdata[2*j] = (temp2[2*j]*cos(phase)+temp2[2*j+1]*sin(phase))*dt/sqrt(2.0*PI);
		Cdata[2*j+1] = (-temp2[2*j]*sin(phase)+temp2[2*j+1]*cos(phase))*dt/sqrt(2.0*PI);
		//Cdata[2*j] = temp2[2*j] * cos(3.0*omega0*chi3/2.0/n0/c*temp*dz/2.0)+temp2[2*j+1] * sin(3.0*omega0*chi3/2.0/n0/c*temp*dz/2.0);
		//Cdata[2*j+1] = -temp2[2*j] * sin(3.0*omega0*chi3/2.0/n0/c*temp*dz/2.0)+temp2[2*j+1] * cos(3.0*omega0*chi3/2.0/n0/c*temp*dz/2.0);
	}

	// FFT
	gsl_fft_complex_forward(Cdata,1,N,wavetable,workspace);

	// temp copy
	for(j=0;j<N;j++)
	{
		temp2[2*j] = Cdata[2*j];
		temp2[2*j+1] = Cdata[2*j+1];
	}

	//delta_z proceed
	for(j=0;j<N;j++)
	{
		omega = start_omega + domega*(double)j - omega0;
		Cdata[2*j] = temp2[2*j] * cos((k2*omega*omega/2.0+k3*omega*omega*omega/6.0)*dz)+temp2[2*j+1] * sin((k2*omega*omega/2.0+k3*omega*omega*omega/6.0)*dz);
		Cdata[2*j+1] = -temp2[2*j] * sin((k2*omega*omega/2.0+k3*omega*omega*omega/6.0)*dz)+temp2[2*j+1] * cos((k2*omega*omega/2.0+k3*omega*omega*omega/6.0)*dz);
	}

	// IFFT
	gsl_fft_complex_inverse(Cdata,1,N,wavetable,workspace);//�t�ϊ�
	
	// temp copy
	for(j=0;j<N;j++)
	{
		temp2[2*j] = Cdata[2*j];
		temp2[2*j+1] = Cdata[2*j+1];
	}

	//z/2 proceed
	for(j=0;j<N;j++)
	{
		temp = temp1[2*j]*temp1[2*j]+temp1[2*j+1]*temp1[2*j+1];
		t = start_t+dt*(double)j;
		phase = 0.5*gamma*temp*dz;
		Cdata[2*j] = temp2[2*j]*cos(phase)+temp2[2*j+1]*sin(phase);
		Cdata[2*j+1] = -temp2[2*j]*sin(phase)+temp2[2*j+1]*cos(phase);
		//Cdata[2*j] = temp2[2*j] * cos(3.0*omega0*chi3/2.0/n0/c*temp*dz/2.0)+temp2[2*j+1] * sin(3.0*omega0*chi3/2.0/n0/c*temp*dz/2.0);
		//Cdata[2*j+1] = -temp2[2*j] * sin(3.0*omega0*chi3/2.0/n0/c*temp*dz/2.0)+temp2[2*j+1] * cos(3.0*omega0*chi3/2.0/n0/c*temp*dz/2.0);
	}

	for(j=0;j<N;j++)
	{
		t = start_t+dt*(double)j;
		phase = (start_omega-omega0)*t;
		temp2[2*j] = (Cdata[2*j]*cos(phase)-Cdata[2*j+1]*sin(phase))*sqrt(2.0*PI)/dt;
		temp2[2*j+1] = (Cdata[2*j]*sin(phase)+Cdata[2*j+1]*cos(phase))*sqrt(2.0*PI)/dt;
	}

	if(i==0)
	{
		sprintf(str,"afterKonishi%d.dat",k);
		fp = fopen(str,"w");
		for(j=0;j<N;j++)
		{
			t = start_t+dt*(double)j;
			fprintf(fp,"%e %e\n",t*1.0e15,c*epsilon0*(temp2[2*j]*temp2[2*j]+temp2[2*j+1]*temp2[2*j+1])/2.0*1.0e-4);
		}
		fclose(fp);
		printf("after %e\n",getEnergy(temp2,dt));
		for(j=0;j<N;j++)
		{
			t = start_t+dt*(double)j;
			phase = (start_omega-omega0)*t;
			temp1[2*j] = (temp2[2*j]*cos(phase)+temp2[2*j+1]*sin(phase))*dt/sqrt(2.0*PI);
			temp1[2*j+1] = (-temp2[2*j]*sin(phase)+temp2[2*j+1]*cos(phase))*dt/sqrt(2.0*PI);
		}
		gsl_fft_complex_forward(temp1,1,N,wavetable,workspace);
		sprintf(str,"afterKonishispectrum%d.dat",k);
		fp = fopen(str,"w");
		for(j=0;j<N;j++)
		{
			omega = start_omega + domega*(double)j;
			fprintf(fp,"%e %e\n",omega*1.0e-15,temp1[2*j]*temp1[2*j]+temp1[2*j+1]*temp1[2*j+1]);
		}
		fclose(fp);
	}

	for(j=0;j<N;j++)
	{
		Cdata[2*j] = temp2[2*j];
		Cdata[2*j+1] = temp2[2*j+1];
	}

}

}//k

gsl_fft_complex_wavetable_free(wavetable);
gsl_fft_complex_workspace_free(workspace);
free(temp1);
free(temp2);
free(Cdata);
return 0;
}

double n_fused_silica(double lambda)
{
	//https://refractiveindex.info
	//Malitson 1965 n 0.21-3.71 um
	double a0,a1,a2;
	double b0,b1,b2;

	a0=0.6961663;
	a1=0.4079426;
	a2=0.8974794;
	b0=0.0684043;
	b1=0.1162414;
	b2=9.896161;

	return sqrt(1.0+a0*lambda*lambda/(lambda*lambda-b0*b0)+a1*lambda*lambda/(lambda*lambda-b1*b1)+a2*lambda*lambda/(lambda*lambda-b2*b2));
}

double n_bk7(double lambda)
{
	//https://refractiveindex.info
	//N-BK7, SCHOTT
	double a0,a1,a2;
	double b0,b1,b2;

	a0=1.03961212;
	a1=0.231792344;
	a2=1.01046945;
	b0=0.00600069867;
	b1=0.0200179144;
	b2=103.560653;

	return sqrt(1.0+a0*lambda*lambda/(lambda*lambda-b0*b0)+a1*lambda*lambda/(lambda*lambda-b1*b1)+a2*lambda*lambda/(lambda*lambda-b2*b2));
}

double getEnergy(double *Cdata, double dt)
{
	int i;
	double energy;
	energy = 0.0;
	for(i=0;i<N;i++)
		 energy += c*epsilon0*(pow(Cdata[2*i],2.0)+pow(Cdata[2*i+1],2.0))/2.0*dt*PI*pow(beam_w,2.0)/2.0;
	return energy;
}

