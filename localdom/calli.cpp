#include "reaction.h"
#include "surface.h"
#include "pot.h"
#include "potPara.h"
#include <fstream>
#include <iostream>
#include "tremendo.h"
#include "structs.h"
#include "definiciones.h"
ofstream misc1("misc1.txt");
ofstream misc2("misc2.txt");
ofstream misc3("misc3.txt");
ofstream misc4("misc4.txt");
ofstream misc5("misc5.txt");
ofstream misc6("misc6.txt");
ofstream misc7("misc7.txt");
ofstream misc8("misc8.txt");
ofstream informe("informe.txt");
using namespace std;
int DOM(int l,double j, double Ecm, double r);
int main(int argc,char* argv[]){
	int n;
	parametros *parm=new struct parametros;
	cout<<"Project managed with Git!!"<<" parameter file: "<<argv[1]<<endl;
       
	if (!parm) Error("No se pudo reservar memoria para parametros");
	const char* input=argv[1];
	LeeParametros(input,parm);
	Capture(parm);
//    if(parm->capture) {Capture(parm);}
    delete[] parm;
	cout<<"Fin! "<<n<<endl;
}
void LeeParametros(const char *fname,struct parametros *x)
{
	char aux[500];
	int potopt,potcm,numst;
	FILE *fp;
	fp = fopen(fname,"r");
	if (!fp) Error("Error al abrir fichero de parametros");
    x->gen12=0;
    x->eikonal=0;
    x->debug=0;
    x->gen_dens_bound=0;
    x->two_trans=0;
    x->one_trans=0;
    x->num_st=-1;
    x->B_numst=-1;
    x->a_numst=-1;
    x->dumb=0;
    x->adiabatico=0;
    x->matching_radio=20.;
    x->relativista=1;
    x->prior=-1;
    x->n_spin=0.5;
    x->remnant=0;
    x->scatt_pot=200;
    x->angle0=0;
    x->angle1=180;
    x->lmin=0;
    x->capture_angular=0;
    strcpy(x->flcoef,"nada");
    strcpy(x->file_dens,"\0");
    potopt=0;
    potcm=0;
    numst=0;
	while(fgets(aux,500,fp)!=NULL)
	{
		if (!potopt) potopt=LeePotencialesOpticos(aux,"InicioPotencialesOpticos",x->pot_opt,fp);
		if (!potcm) potcm=LeePotencialesCampoMedio(aux,"InicioCampoMedio",x->pot,fp);
		if (!numst && (x->two_trans || x->knockout || x->one_trans || x->capture)) numst=LeeEstados(aux,"InicioEstados",x->st,fp);
		ReadParS(aux,"flcoef",x->flcoef);
		ReadParS(aux,"fl_log",x->fl_log);
		ReadParS(aux,"file_dens",x->file_dens);
		ReadParS(aux,"fl_energias",x->fl_energias);
		ReadParS(aux,"fl_formfactor",x->fl_formfactor);
		ReadParS(aux,"fl_funondas",x->fl_funondas);
		ReadParS(aux,"fl_potcm",x->fl_potcm);
		ReadParS(aux,"fl_amplitudes",x->fl_amplitudes);
		ReadParS(aux,"fl_cross_succ",x->fl_cross_succ);
		ReadParS(aux,"fl_cross_sim",x->fl_cross_sim);
		ReadParS(aux,"fl_cross_non",x->fl_cross_non);
		ReadParS(aux,"fl_cross_tot",x->fl_cross_tot);
		ReadParS(aux,"fl_fundamental",x->fl_fundamental);
		ReadParS(aux,"fl_diagrama",x->fl_diagrama);
		ReadParS(aux,"fl_espectro",x->fl_espectro);
		ReadParS(aux,"fl_dw",x->fl_dw);
		ReadParS(aux,"fl_gf",x->fl_gf);
		ReadParS(aux,"a_tipo_fun",x->a_tipo_fun);
		ReadParS(aux,"B_tipo_fun",x->B_tipo_fun);
		ReadParD(aux,"a_potcm",&(x->a_potcm));
                ReadParF(aux,"enerange_min",&(x->enerange_min));
                ReadParF(aux,"enerange_max",&(x->enerange_max));
                ReadParF(aux,"enerange_step",&(x->enerange_step));
		ReadParS(aux,"proyectil",x->proyectil);
		ReadParS(aux,"unidades",x->unidades);
		ReadParD(aux,"num_st",&(x->num_st));
		ReadParD(aux,"dompot_n",&(x->dompot_n));
		ReadParD(aux,"dompot_p",&(x->dompot_p));
		ReadParD(aux,"capture_angular",&(x->capture_angular));
		ReadParD(aux,"two_trans",&(x->two_trans));
		ReadParD(aux,"one_trans",&(x->one_trans));
		ReadParD(aux,"capture",&(x->capture));
		ReadParD(aux,"knockout",&(x->knockout));
		ReadParD(aux,"zerorange",&(x->zerorange));
		ReadParD(aux,"relativista",&(x->relativista));
		ReadParD(aux,"eikonal",&(x->eikonal));
		ReadParD(aux,"twonuceikonal",&(x->twonuceikonal));
		ReadParD(aux,"optico_strip",&(x->optico_strip));
		ReadParD(aux,"optico_dif",&(x->optico_dif));
		ReadParD(aux,"b_puntos",&(x->b_puntos));
		ReadParD(aux,"k_t_puntos",&(x->k_t_puntos));
		ReadParD(aux,"remnant",&(x->remnant));
		ReadParD(aux,"core_pot",&(x->core_pot));
		ReadParD(aux,"scatt_pot",&(x->scatt_pot));
		ReadParF(aux,"emin",&(x->emin));
		ReadParF(aux,"n_spin",&(x->n_spin));
		ReadParF(aux,"emax",&(x->emax));
		ReadParD(aux,"puntos",&(x->puntos));
		ReadParD(aux,"adiabatico",&(x->adiabatico));
		ReadParD(aux,"prior",&(x->prior));
		ReadParF(aux,"radio",&(x->radio));
		ReadParF(aux,"matching_radio",&(x->matching_radio));
		ReadParF(aux,"delta_r",&(x->delta_r));
		ReadParD(aux,"gen12",&(x->gen12));
		ReadParD(aux,"gen_dens_bound",&(x->gen_dens_bound));
		ReadParD(aux,"debug",&(x->debug));
		ReadParD(aux,"lmin",&(x->lmin));
		ReadParD(aux,"lmax",&(x->lmax));
		ReadParD(aux,"ltransfer",&(x->ltransfer));
		ReadParD(aux,"num_cm",&(x->num_cm));
		ReadParD(aux,"num_opt",&(x->num_opt));
		ReadParD(aux,"id_pot_dens",&(x->id_pot_dens));
		ReadParD(aux,"B_numst",&(x->B_numst));
		ReadParD(aux,"a_numst",&(x->a_numst));
		ReadParD(aux,"B_potcm",&(x->B_potcm));
		ReadParD(aux,"base1",&(x->base1));
		ReadParD(aux,"base2",&(x->base2));
		ReadParD(aux,"base3",&(x->base3));
		ReadParD(aux,"dumb",&(x->dumb));
		ReadParD(aux,"lmax_PT",&(x->lmax_PT));
		ReadParD(aux,"lmax_RN",&(x->lmax_RN));
		ReadParF(aux,"m_A",&(x->m_A));
		ReadParF(aux,"m_a",&(x->m_a));
		ReadParF(aux,"m_B",&(x->m_B));
		ReadParF(aux,"m_b",&(x->m_b));
		ReadParF(aux,"Z_A",&(x->Z_A));
		ReadParF(aux,"Z_a",&(x->Z_a));
		ReadParF(aux,"Z_B",&(x->Z_A));
		ReadParF(aux,"Z_b",&(x->Z_a));
		ReadParF(aux,"energia_lab",&(x->energia_lab));
		ReadParF(aux,"Qvalue",&(x->Qvalue));
		ReadParF(aux,"int_Qvalue",&(x->int_Qvalue));
		ReadParF(aux,"J_a",&(x->J_a));
		ReadParF(aux,"J_A",&(x->J_A));
		ReadParF(aux,"J_b",&(x->J_b));
		ReadParF(aux,"J_B",&(x->J_B));
		ReadParF(aux,"dw_spinA",&(x->dw_spinA));
		ReadParF(aux,"dw_spinB",&(x->dw_spinB));
		ReadParF(aux,"P_masa",&(x->P_masa));
		ReadParF(aux,"T_masa",&(x->T_masa));
		ReadParF(aux,"n1_masa",&(x->n1_masa));
		ReadParF(aux,"n2_masa",&(x->n2_masa));
		ReadParF(aux,"res_masa",&(x->res_masa));
		ReadParF(aux,"P_carga",&(x->P_carga));
		ReadParF(aux,"T_carga",&(x->T_carga));
		ReadParF(aux,"P_N",&(x->P_N));
		ReadParF(aux,"T_N",&(x->T_N));
		ReadParF(aux,"res_N",&(x->res_N));
		ReadParF(aux,"n1_carga",&(x->n1_carga));
		ReadParF(aux,"n2_carga",&(x->n2_carga));
		ReadParF(aux,"n1_N",&(x->n1_N));
		ReadParF(aux,"res_carga",&(x->res_carga));
		ReadParD(aux,"folding",&(x->folding));
		ReadParD(aux,"fermi_dens",&(x->fermi_dens));
		ReadParD(aux,"gauss_dens",&(x->gauss_dens));
		ReadParD(aux,"koning_delaroche",&(x->koning_delaroche));
		ReadParF(aux,"P_sigma_dens",&(x->P_sigma_dens));
		ReadParF(aux,"T_sigma_dens",&(x->T_sigma_dens));
		ReadParF(aux,"P_radio_dens",&(x->P_radio_dens));
		ReadParF(aux,"T_radio_dens",&(x->T_radio_dens));
		ReadParF(aux,"lambda",&(x->lambda));
		ReadParF(aux,"angle0",&(x->angle0));
		ReadParF(aux,"angle1",&(x->angle1));
		ReadParD(aux,"pot_transfer",&(x->pot_transfer));
		ReadParD(aux,"optico_ingreso",&(x->optico_ingreso));
		ReadParD(aux,"optico_intermedio",&(x->optico_intermedio));
		ReadParD(aux,"optico_salida",&(x->optico_salida));
		ReadParD(aux,"successive",&(x->successive));
		ReadParD(aux,"simultaneous",&(x->simultaneous));
		if(x->B_numst>0) ReadParMultD(aux,"B_estados",x->B_numst,x->B_estados);
		if(x->a_numst>0) ReadParMultD(aux,"a_estados",x->a_numst,x->a_estados);
		ReadParF(aux,"V0pairing",&(x->V0pairing));
		ReadParF(aux,"Vrpairing",&(x->Vrpairing));
		ReadParF(aux,"rho0",&(x->rho0));
		ReadParF(aux,"pexp",&(x->pexp));
		ReadParF(aux,"r_Ccmin",&(x->r_Ccmin));
		ReadParF(aux,"r_Ccmax",&(x->r_Ccmax));
		ReadParF(aux,"r_A2min",&(x->r_A2min));
		ReadParF(aux,"r_A2max",&(x->r_A2max));
		ReadParD(aux,"rCc_puntos",&(x->rCc_puntos));
		ReadParD(aux,"rA2_puntos",&(x->rA2_puntos));
		ReadParD(aux,"theta_puntos",&(x->theta_puntos));
		ReadParD(aux,"cross_puntos",&(x->cross_puntos));
		ReadParD(aux,"form_factor",&(x->form_factor));
		ReadParF(aux,"a_Sn",&(x->a_Sn));
		ReadParF(aux,"B_Sn",&(x->B_Sn));
	}
//	cout<<"x->a_estados: "<<x->a_estados[0]<<endl;
//	exit(0);
	if((x->a_numst+x->B_numst)>(x->num_st)) Error("El numero de estados total debe ser la suma de los estados de ambos n�cleos");
	if(potopt!=x->num_opt) Error("El numero de potenciales opticos leidos no coincide con el especificado");
	if(potcm!=x->num_cm) Error("El numero de potenciales de campo medio leidos no coincide con el especificado");
	if(numst!=x->num_st && x->two_trans && !(!strcmp(x->a_tipo_fun,"li")) && !(!strcmp(x->B_tipo_fun,"li")))
		Error("El numero de estados leidos no coincide con el especificado");
	if(x->prior==-1 && x->two_trans==1) Error("La representacion (post-prior) no ha sido definida");
	fclose(fp);
}
void Capture(struct parametros* parm)
{
	cout<<"********************************************************************************"<<endl;
	cout<<"*                                                                              *"<<endl;
	cout<<"*                       CAPTURA NEUTRONICA                                     *"<<endl;
	cout<<"*                                                                              *"<<endl;
	cout<<"********************************************************************************"<<endl;
	int n,m,indx_pot_a,indx_pot_B,indx_st,indx_ingreso,indx_intermedio,indx_salida,indx_core,indx_transfer,indx_scatt;
	double energia,etrial,vmax,vmin,energia_ws,absorcion,carga_out,carga_trans,A_P,A_T;
	cout<<" Energia laboratorio **********************"<<parm->energia_lab<<endl;
	InicializaOneTrans(parm);
	double* D0=new double[1];
	double* rms=new double[1];
	cout<<"Generando potenciales de campo medio"<<endl;
	HanShiShen(parm->energia_lab,parm->T_N,parm->T_carga);
	for(n=0;n<parm->num_cm;n++)
	{
		GeneraPotencialCM(parm,&(parm->pot[n]));
		if(parm->a_potcm==parm->pot[n].id) indx_pot_a=n;
		if(parm->B_potcm==parm->pot[n].id) indx_pot_B=n;
	}
	for (n=0;n<parm->num_opt;n++)
	{
		if(parm->optico_ingreso==parm->pot_opt[n].id) indx_ingreso=n;
		if(parm->optico_salida==parm->pot_opt[n].id) indx_salida=n;
		if(parm->optico_intermedio==parm->pot_opt[n].id) indx_intermedio=n;
		if(parm->core_pot==parm->pot_opt[n].id) indx_core=n;
		if(parm->pot_transfer==parm->pot_opt[n].id) indx_transfer=n;
	}
	A_P=parm->P_N+parm->P_carga;
	A_T=parm->T_N+parm->T_carga;
	cout<<"Generando potenciales opticos"<<endl;
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_ingreso]),A_T,A_P);
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_salida]),parm->m_B,parm->m_b);
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_intermedio]),parm->m_a-1.,parm->m_A);
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_core]),parm->m_a-1.,parm->m_A);
	GeneraPotencialOptico(parm,&(parm->pot_opt[indx_transfer]),1.,parm->m_A);


	cout<<"Masa del nucleo compuesto: "<<parm->res_masa<<endl;
	carga_trans=parm->res_carga-parm->T_carga;
	cout<<"Carga de la particula transferida: "<<carga_trans<<endl;
	carga_out=parm->P_carga-carga_trans;
	cout<<"Carga de la particula emitida: "<<carga_out<<endl;
	cout<<"Masa  de la particula transferida: "<<parm->n1_masa<<endl;
	cout<<"Masa  de la particula emitida: "<<parm->m_b<<endl;
	cout<<"Generando el estado del nucleo a"<<endl;
	/* Genera niveles del n�cleo 'a' */
	for (n=0;n<parm->a_numst;n++)
	{
		for(m=0;m<parm->num_st;m++)
		{
			if(parm->a_estados[n]==parm->st[m].id) indx_st=m;
		}
		if(parm->st[indx_st].energia<0.)
			GeneraEstadosPI(&(parm->pot[indx_pot_a]),&(parm->st[indx_st]),parm->radio,parm->puntos,carga_trans*(carga_out),parm,1,
					parm->n1_masa*parm->m_b/(parm->n1_masa+parm->m_b),D0,rms);
//			HulthenWf(&(parm->st[indx_st]),parm->radio,parm->puntos);
		else
		{
			GeneraEstadosContinuo(&(parm->pot_opt[indx_scatt]),&(parm->st[indx_st]),parm->radio,parm->puntos,carga_trans*(carga_out)
					,parm,parm->n1_masa*parm->m_b/(parm->n1_masa+parm->m_b));
		}
	}
	cout<<"Generando niveles nucleo B"<<endl;
	/* Genera niveles del n�cleo 'B' */
	for (n=0;n<parm->B_numst;n++)
	{
		for(m=0;m<parm->num_st;m++)
		{
			if(parm->B_estados[n]==parm->st[m].id) indx_st=m;
		}
		cout<<"Masa reducida: "<<parm->m_B*parm->n1_masa/(parm->m_B+parm->n1_masa)<<endl;
		if(parm->st[indx_st].energia<0.)
			GeneraEstadosPI(&(parm->pot[indx_pot_B]),&(parm->st[indx_st]),parm->radio,parm->puntos,
					carga_trans*parm->T_carga,parm,0,parm->m_B*parm->n1_masa/(parm->m_B+parm->n1_masa),D0,rms);
		else
		{
			GeneraEstadosContinuo(&(parm->pot_opt[indx_scatt]),&(parm->st[indx_st]),parm->radio,parm->puntos,
					carga_trans*parm->T_carga,parm,parm->m_B*parm->n1_masa/(parm->m_B+parm->n1_masa));
		}
		absorcion=Absorcion2(&(parm->pot_opt[indx_intermedio]),&(parm->st[indx_st]));
	}
	cout<<"Profundidad pozo: "<<parm->pot[indx_pot_B].V<<endl;
	EscribeEstados(parm->puntos,parm->st,parm->num_st,parm);
	EscribePotencial(parm->puntos,parm->pot,parm->num_cm,parm);
	EscribePotencialOptico(parm->puntos,parm->pot_opt,parm->num_opt,parm);
//	cout<<"koning: "<<parm->koning_delaroche<<endl;
	if(parm->koning_delaroche==2) AmplitudeCaptureCC(parm);
	else AmplitudeCapture(parm);
}
void AmplitudeCapture(struct parametros* parm)
{

	parametros_integral *dim1=new parametros_integral;
	parametros_integral *dim2=new parametros_integral;
	parametros_integral *dim3=new parametros_integral;
	parametros_integral *dim4=new parametros_integral;
	complejo* exp_delta_coulomb_i=new complejo[parm->lmax];
	complejo* exp_delta_coulomb_f=new complejo[parm->lmax];
	double eta_f=parm->Z_a*parm->Z_A*E2HC*parm->mu_Bb*AMU/(HC*parm->k_Bb);
	double eta_i=parm->eta;
	double step,rn,energia_out,energia_trans,k_p,k_n,cross,elastic_cross,
	theta,costheta,D0,rhoE,sigma_const,escala,r_source,velocidad,
	cross_total,cross_total_elasticb,redfac,r_F,absorcion,e_res,rhoE_n,N_A,
	carga_out,carga_trans,km,rAn,Ecm,Ecm_out,cross_total_breakup;
	distorted_wave* fl=new distorted_wave;
	distorted_wave* gl_up=new distorted_wave;
	distorted_wave* gl_down=new distorted_wave;
	
	distorted_wave* funcion_regular_up=new distorted_wave[2];
	distorted_wave* funcion_irregular_up=new distorted_wave[2];
	distorted_wave* funcion_regular_down=new distorted_wave[2];
	distorted_wave* funcion_irregular_down=new distorted_wave[2];
	potencial_optico *optico=new potencial_optico;
	potencial_optico *core=new potencial_optico;
	potencial_optico* v_up=new potencial_optico[1];
	potencial_optico* v_down=new potencial_optico[1];
	potencial_optico* vp_up=new potencial_optico[1];
	potencial_optico* vp_down=new potencial_optico[1];
	potencial_optico* pot_dumb=new potencial_optico;
	estado* st=new estado;
	estado* st_fin=new estado;
	ofstream fp1("dw_out1trans.txt");
	ofstream fp2("dw_in1trans.txt");
	ofstream fp3;
	fp3.open("talys1.txt");
	ofstream fp4;
	fp4.open("talys2.txt");
	ofstream fp5;
	fp5.open("talys_angular1.txt");
	ofstream fp6;
	fp6.open("talys_angular2.txt");
	ofstream fp7;
	string name40("40Ca_");
	string name48("48Ca_");
	string name60("60Ca_");
	string name40Spin("40Ca_");
	string name48Spin("48Ca_");
	string name60Spin("60Ca_");
	string name40Angular("40Ca_");
	string name48Angular("48Ca_");
	string name60Angular("60Ca_");
	stringstream ss;
	ss<<parm->energia_lab;
	string name_energy=ss.str();
	//char* name_energy=itoa(int(ceil(parm->energia_lab)));
	//itoa(parm->energia_lab,name_energy,10);
	name40+=name_energy+"MeV";
	name48+=name_energy+"MeV";
	name60+=name_energy+"MeV";
	name40Spin+=name_energy+"MeVSpinParity";
	name48Spin+=name_energy+"MeVSpinParity";
	name60Spin+=name_energy+"MeVSpinParity";
	name40Angular+=name_energy+"MeVAngular";
	name48Angular+=name_energy+"MeVAngular";
	name60Angular+=name_energy+"MeVAngular";
	const char * n40 = name40.c_str();
	const char * n48 = name48.c_str();
	const char * n60 = name60.c_str();

	const char * n40S = name40Spin.c_str();
	const char * n48S = name48Spin.c_str();
	const char * n60S = name60Spin.c_str();

	const char * n40A = name40Angular.c_str();
	const char * n48A = name48Angular.c_str();
	const char * n60A = name60Angular.c_str();

	if (parm->dompot_n==0) fp7.open(n40S);
	if (parm->dompot_n==2) fp7.open(n48S);
	if (parm->dompot_n==4) fp7.open(n60S);
	if (parm->dompot_n>49) fp7.open("SpinParity.txt");
	ofstream fp8;
	fp8.open("Jutta_angular.txt");
	ofstream fp9;
	if (parm->dompot_n==0) fp9.open(n40);
	if (parm->dompot_n==2) fp9.open(n48);
	if (parm->dompot_n==4) fp9.open(n60);
	if (parm->dompot_n>49) fp9.open("dsdE.txt");
	ofstream fp10;
	if (parm->dompot_n==0) fp10.open(n40A);
	if (parm->dompot_n==2) fp10.open(n48A);
	if (parm->dompot_n==4) fp10.open(n60A);
	if (parm->dompot_n>49) fp10.open("dsdEdO.txt");
	if (parm->dompot_n>49) {parm->dompot_n=0;parm->dompot_p=1;}
	complejo* S=new complejo[parm->lmax];
	complejo**** rho=tensor4_cmpx(parm->rCc_puntos,parm->lmax,parm->lmax+1,parm->lmax);
	complejo* rho_test=new complejo [parm->puntos];
	complejo* rhom=new complejo[parm->lmax+1];
	complejo**** non=tensor4_cmpx(parm->rCc_puntos,parm->lmax,parm->lmax+1,parm->lmax);
	complejo**** dumb=tensor4_cmpx(parm->rCc_puntos,parm->lmax,parm->lmax+1,parm->lmax);
	complejo* non_test=new complejo [parm->puntos];
	complejo* nonm=new complejo[parm->lmax+1];
	complejo**** phi_up=tensor4_cmpx(parm->rCc_puntos,parm->lmax,parm->lmax+1,parm->lmax);
	complejo**** phi_down=tensor4_cmpx(parm->rCc_puntos,parm->lmax,parm->lmax+1,parm->lmax);
	complejo**** phi_resonant=tensor4_cmpx(parm->rCc_puntos,parm->lmax,parm->lmax+1,parm->lmax);
	complejo*** Teb=tensor_cmpx(parm->lmax,parm->lmax+1,parm->lmax);
	complejo** GreenFunction=matriz_cmpx(300,300);
	complejo* phi_test=new complejo [parm->puntos];
	complejo* phi_res=new complejo [parm->puntos];
	complejo* phim=new complejo[parm->lmax+1];
	complejo* wf=new complejo[1000];
	double* rg=new double[300];
	complejo pot_p;
	complejo pot_n;
    pot opt_in(parm->dompot_n);
    pot opt_out(parm->dompot_p);
	double* total_break=new double[parm->lmax+1];
	double* inc_break=new double[parm->lmax+1];
	double* inc_break_lmas=new double[parm->lmax+1];
	double* inc_break_lmenos=new double[parm->lmax+1];
	double* cross_up=new double[parm->lmax+1];
	double* cross_down=new double[parm->lmax+1];
	double* elastic_break=new double[parm->lmax+1];
	double* direct=new double[1];
	double* non_orth=new double[1];
	double* cross_term=new double[1];
	double** Al=matriz_dbl(2*parm->lmax+1,parm->lmax);
	int l,lp,ld,indx_salida,indx_ingreso,indx_core,indx_neutron_target,indx_st,n,la,m,len,flag,n1;
	complejo rhofac,ampli,wronskiano,wronskiano_up,wronskiano_down,fl_int,gl_int,fl_source,gl_source,st_source,st_int,lorentz,is_pot_im,is_pot_im_out;
	dim1->a=parm->r_Ccmin;
	dim1->b=parm->r_Ccmax;
	dim1->num_puntos=parm->rCc_puntos;
	dim2->a=0.;
	dim2->b=PI;
	dim2->num_puntos=parm->theta_puntos;
	dim3->a=parm->r_A2min;
	dim3->b=parm->r_A2max;
	dim3->num_puntos=parm->rA2_puntos;
	dim4->a=parm->r_A2min;
	dim4->b=parm->r_A2max;
	dim4->num_puntos=parm->rA2_puntos;
	GaussLegendre(dim1->puntos,dim1->pesos,dim1->num_puntos);
	GaussLegendre(dim2->puntos,dim2->pesos,dim2->num_puntos);
	GaussLegendre(dim3->puntos,dim3->pesos,dim3->num_puntos);
	GaussLegendre(dim4->puntos,dim4->pesos,dim4->num_puntos);
	D0=10.;
    redfac=2.*AMU/(HC*HC);
	cout<<"Masa del proyectil: "<<parm->P_masa<<endl;
	cout<<"Masa del blanco: "<<parm->T_masa<<endl;
	cout<<"Masa del nucleo compuesto: "<<parm->res_masa<<endl;
	carga_trans=parm->res_carga-parm->T_carga;
	cout<<"Carga de la particula transferida: "<<carga_trans<<endl;
	carga_out=parm->P_carga-carga_trans;
	cout<<"Carga de la particula emitida: "<<carga_out<<endl;
	cout<<"Masa  de la particula transferida: "<<parm->n1_masa<<endl;
	cout<<"Masa reducida de la particula emitida: "<<parm->m_b<<endl;
	km=(parm->m_A+1.)/parm->m_A;

	/*Selecciona los potenciales opticos en los distintos canales*/
	for (n=0;n<parm->num_opt;n++)
	{
		if(parm->optico_ingreso==parm->pot_opt[n].id) indx_ingreso=n;
		if(parm->optico_salida==parm->pot_opt[n].id) indx_salida=n;
		if(parm->optico_intermedio==parm->pot_opt[n].id) indx_neutron_target=n;
		if(parm->core_pot==parm->pot_opt[n].id) indx_core=n;
		//if(parm->pot_transfer==parm->pot_opt[n].id) v_up=&(parm->pot_opt[n]);
		//	if(parm->scatt_pot==parm->pot_opt[n].id) v_down=&(parm->pot_opt[n]);
	}
	//GeneraPotencialOptico(parm,v_up,1.,parm->m_A);
	//GeneraPotencialOptico(parm,v_down,1.,parm->m_A);
	v_up->puntos=parm->puntos;
	v_down->puntos=parm->puntos;
	vp_up->puntos=parm->puntos;
	vp_down->puntos=parm->puntos;	
	v_up->radio=parm->radio;
	v_down->radio=parm->radio;
	vp_up->radio=parm->radio;
	vp_down->radio=parm->radio;
	v_up->radio_coul=5.;
	v_down->radio_coul=5.;
	vp_up->radio_coul=5.;
	vp_down->radio_coul=5.;
	v_up->Vso=5.;
	v_down->Vso=5.;
	vp_up->Vso=5.;
	vp_down->Vso=5.;
	v_up->radioso=5.;
	v_down->radioso=5.;
	vp_up->radioso=5.;
	vp_down->radioso=5.;
	v_up->aso=0.65;
	v_down->aso=0.65;
	vp_up->aso=0.65;
	vp_down->aso=0.65;
	cout<<"Energia de centro de masa: "<<parm->energia_cm<<endl;
	cout<<"Momento inicial: "<<parm->k_Aa<<endl;
	cout<<"Momento final: "<<parm->k_Bb<<endl;
	cout<<"Masa reducida deuteron: "<<parm->mu_Aa<<endl;
	cout<<"Masa reducida proton: "<<parm->mu_Bb<<endl;
	cout<<"Masa reducida neutron: "<<parm->m_A/(parm->m_A+1.)<<endl;
	/*Calculo de las amplitudes de transferencia**************************************************************************/
	for(n=0;n<parm->num_st;n++)
	{
		if (parm->a_estados[0] == parm->st[n].id) st= &(parm->st[n]);
		if (parm->B_estados[0] == parm->st[n].id) st_fin= &(parm->st[n]);
	}
	step=double(parm->radio/parm->puntos);
	absorcion=Absorcion2(&(parm->pot_opt[indx_neutron_target]),st_fin);
	cout<<"Absorcion: "<<absorcion<<endl;
	for(la=0;la<parm->lmax;la++)
	{
		exp_delta_coulomb_i[la]=exp(I*(deltac(la,eta_i)));
		exp_delta_coulomb_f[la]=exp(I*(deltac(la,eta_f)));
	}
	velocidad=C*sqrt(2*parm->energia_lab/(2.*AMU));
	sigma_const=2.*parm->mu_Aa*AMU/(HC*HC*parm->k_Aa);
	len=strlen(parm->unidades);
	if(!strncmp(parm->unidades,"milib",len)) flag=1;
	if(!strncmp(parm->unidades,"fm2",len)) flag=2;
	if(!strncmp(parm->unidades,"b",len)) flag=3;
	if(!strncmp(parm->unidades,"microb",len)) flag=4;
	switch(flag)
	{
	case 1:
		escala=10.;
		cout<<"Seccion efi caz medida en milibarn"<<endl;
		break;
	case 2:
		escala=1.;
		cout<<"Seccion eficaz medida en fm^2"<<endl;
		break;
	case 3:
		escala=0.01;
		cout<<"Seccion eficaz medida en barn"<<endl;
		break;
	case 4:
		escala=10000.;
		cout<<"Seccion eficaz medida en microbarn"<<endl;
		break;
	default:
		Error("Unidades desconocidas para la secci�n eficaz");
		break;
	}

	if(parm->koning_delaroche==1) cout<<"*****************************************************"<<endl<<
									    "***** Potencial neutron-blanco Koning-Delaroche *****"<<endl<<
									    "*****************************************************"<<endl;
	if(parm->koning_delaroche==0) cout<<"*****************************************************"<<endl<<
									    "***** Using DOM potentials *****"<<endl<<
									    "*****************************************************"<<endl;
        if(parm->koning_delaroche<0 || parm->koning_delaroche>1) {cout<<"Error!: koning_delaroche has to be either 0 or 1  "<<endl; exit(0);}
	r_F=1000.;
	e_res=st_fin->energia;
	for(energia_out=parm->enerange_min;energia_out<parm->enerange_max;energia_out+=parm->enerange_step)
	{
		Ecm_out=((parm->T_masa)*energia_out/(parm->n1_masa+(parm->T_masa)));
		Ecm=parm->energia_cm-Ecm_out-2.2245;
		energia_trans=(parm->n1_masa+parm->T_masa)*Ecm/(parm->T_masa);
		cout<<endl<<endl<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
		cout<<"Energy of detected cluster: "<<energia_out<<"  "<<"Energy of absorbed cluster: "
				<<energia_trans<<endl;
		fp9<<energia_out<<"  "<<energia_trans<<"  ";
		misc1<<"& Energy of detected cluster: "<<energia_out<<"    Energy of absorbed cluster: "<<energia_trans<<endl;
//		misc2<<endl<<"*********************  Ep= "<<energia_out<<" ****************************"<<endl;
		k_n=sqrt(2.*parm->n1_masa*AMU*Ecm)/HC;
		k_p=sqrt(2.*parm->m_b*AMU*Ecm_out)/HC;
		rhoE=parm->m_b*AMU*k_p/(8.*PI*PI*PI*HC*HC);
		rhoE_n=parm->n1_masa*AMU*k_n/(8.*PI*PI*PI*HC*HC);
		eta_f=carga_out*parm->res_carga*E2HC*(parm->m_b*parm->T_masa/(parm->m_b+parm->T_masa))*AMU/(HC*k_p);
		cross_total=0.;
		cross_total_elasticb=0.;
		opt_in.potential(1.,0,0,energia_trans);
		opt_out.potential(1.,0,0,energia_out);
		opt_in.GetPot(0.5,0,0.5);
		opt_out.GetPot(0.5,0,0.5);
		is_pot_im_out=opt_out.Imag;
		is_pot_im=opt_in.Imag;
		for(l=parm->lmin;l<parm->ltransfer;l++)
		{
		  if(parm->koning_delaroche==0){
		 for(n=0;n<parm->puntos;n++)  
		   {
  			rn=step*(n+1.);
  			opt_in.GetPot(rn,l,l+0.5);
  			v_up->r[n]=rn;
  			v_up->pot[n]=opt_in.Real+I*opt_in.Imag;
  			if(abs(is_pot_im)<1.e-10)  v_up->pot[n]+=I*0.001*opt_in.Real;			      
  			if(l>0) opt_in.GetPot(rn,l,l-0.5);
  			v_down->r[n]=rn;
  			v_down->pot[n]=opt_in.Real+I*opt_in.Imag;
  			if(abs(is_pot_im)<1.e-10) v_down->pot[n]+=I*0.001*opt_in.Real;
		    }
		  }
		if(parm->koning_delaroche==1){
				for(n=0;n<parm->puntos;n++)
				{
					rn=step*(n+1.);
					KoningDelaroche(energia_out,parm->T_N,parm->T_carga,rn,&pot_p,
							&pot_n,l,l+0.5,vp_down,pot_dumb);
					vp_up->r[n]=rn;
					vp_up->pot[n]=pot_p;
					if(l>0) KoningDelaroche(energia_out,parm->T_N,parm->T_carga,rn,&pot_p,
							&pot_n,l,l-0.5,vp_down,pot_dumb);
					if(l==0) KoningDelaroche(energia_out,parm->T_N,parm->T_carga,rn,&pot_p,
							&pot_n,l,0.5,vp_down,pot_dumb);
					vp_down->r[n]=rn;
					vp_down->pot[n]=pot_p;
					KoningDelaroche(energia_trans,parm->T_N,parm->T_carga,rn,&pot_p,
							&pot_n,l,l+0.5,pot_dumb,v_down);
					v_up->r[n]=rn;
					v_up->pot[n]=pot_n;
					v_down->pot[n]=pot_n;
					if(l>0) KoningDelaroche(energia_trans,parm->T_N,parm->T_carga,rn,&pot_p,
							&pot_n,l,l-0.5,pot_dumb,v_down);
					if(l==0) KoningDelaroche(energia_out,parm->T_N,parm->T_carga,rn,&pot_p,
												&pot_n,l,0.5,pot_dumb,v_down);
					v_down->r[n]=rn;
					v_down->pot[n]=pot_n;
				}
			}
			cout<<"L: "<<l<<endl;
			funcion_regular_up[0].energia=Ecm;
			funcion_irregular_up[0].energia=Ecm;
			funcion_regular_up[0].l=l;
			funcion_irregular_up[0].l=l;
			funcion_regular_up[0].j=l+parm->n_spin;     
			funcion_irregular_up[0].j=l+parm->n_spin;
			funcion_regular_up[1].energia=Ecm;
			funcion_irregular_up[1].energia=Ecm;
			funcion_regular_up[1].l=l;
			funcion_irregular_up[1].l=l;
			funcion_regular_up[1].j=l+parm->n_spin;
			funcion_irregular_up[1].j=l+parm->n_spin;
			if (Ecm<=0.) wronskiano_up=GeneraGreenFunctionLigada(&(funcion_regular_up[0]),&(funcion_irregular_up[0]),
					v_up,parm->radio,parm->puntos,0.,
					parm->n1_masa*parm->T_masa/(parm->n1_masa+parm->T_masa),parm->n_spin);
			if (Ecm>0.)
			{
				GeneraGreenFunction(funcion_regular_up,funcion_irregular_up,v_up,
						0.,parm->n1_masa*parm->T_masa/(parm->n1_masa+parm->T_masa),
						parm->radio,parm->puntos,parm->matching_radio,parm->n_spin);
				wronskiano_up=k_n;
			}
			funcion_regular_down[1].energia=Ecm;
			funcion_irregular_down[1].energia=Ecm;
			funcion_regular_down[1].l=l;
			funcion_irregular_down[1].l=l;
			funcion_regular_down[1].j=l-parm->n_spin;
			funcion_irregular_down[1].j=l-parm->n_spin;
			funcion_regular_down[0].energia=Ecm;
			funcion_irregular_down[0].energia=Ecm;
			funcion_regular_down[0].l=l;
			funcion_irregular_down[0].l=l;
			funcion_regular_down[0].j=l-parm->n_spin;
			funcion_irregular_down[0].j=l-parm->n_spin;
			if(l==0) {
				funcion_irregular_down[0].j=parm->n_spin; funcion_regular_down[0].j=parm->n_spin;
				funcion_irregular_down[1].j=parm->n_spin; funcion_regular_down[1].j=parm->n_spin;
			}
			if (Ecm<=0.) wronskiano_down=GeneraGreenFunctionLigada
					(&(funcion_regular_down[1]),&(funcion_irregular_down[1]),
							v_down,parm->radio,parm->puntos,0.,
							parm->n1_masa*parm->T_masa/(parm->n1_masa+parm->T_masa),parm->n_spin);
			if (Ecm>0.)
			{
				GeneraGreenFunction(funcion_regular_down,funcion_irregular_down,v_down,
						0.,parm->n1_masa*parm->T_masa/(parm->n1_masa+parm->T_masa),
						parm->radio,parm->puntos,parm->matching_radio,parm->n_spin);
				wronskiano_down=k_n;
			}
			for(lp=0;lp<parm->lmax;lp++)
			{
				if(parm->koning_delaroche==0){
					for(n=0;n<parm->puntos;n++)
					{
						rn=step*(n+1.);
						opt_out.GetPot(rn,lp,lp+0.5);
						vp_up->r[n]=rn;
						vp_up->pot[n]=opt_out.Real+I*opt_out.Imag;
						if(abs(is_pot_im_out)<1.e-10)  vp_up->pot[n]+=I*0.001*opt_out.Real;
						if(lp>0) opt_out.GetPot(rn,lp,lp-0.5);
						vp_down->r[n]=rn;
						vp_down->pot[n]=opt_out.Real+I*opt_out.Imag;
						if(abs(is_pot_im_out)<1.e-10) vp_down->pot[n]+=I*0.005*opt_out.Real;
					}
				}
			if(parm->remnant==1 && parm->prior==1) {
			  GeneraRemnant(optico,core,&parm->pot_opt[indx_ingreso],vp_down,parm->T_carga*parm->P_carga,
							0.,0,0,parm->mu_Aa,parm->m_b);

			}
				gl_up->energia=Ecm_out;
				gl_up->l=lp;
				gl_up->spin=parm->n_spin;
				gl_up->j=lp+parm->n_spin;
				GeneraDWspin(gl_up,vp_up,0.,parm->m_b*parm->res_masa/(parm->m_b+parm->res_masa),
					     parm->radio,parm->puntos,parm->matching_radio,&fp2);
				gl_down->energia=Ecm_out;
				gl_down->l=lp;
				gl_down->spin=parm->n_spin;
				gl_down->j=lp-parm->n_spin;
				if(lp==0) gl_down->j=lp;
				GeneraDWspin(gl_down,vp_down,0.,parm->m_b*parm->res_masa/(parm->m_b+parm->res_masa),
						parm->radio,parm->puntos,parm->matching_radio,&fp2);
				for(n=0;n<dim1->num_puntos;n++){
					for(m=0;m<=lp;m++){
						rho[n][l][m][lp]=0.;
						non[n][l][m][lp]=0.;
					}
				}
				exp_delta_coulomb_f[lp]=exp(I*(deltac(lp,eta_f)));
				for(ld=abs(l-lp);(ld<=l+lp)&&(ld<parm->lmax);ld++)
				{
					rhofac=(16.*pow(PI,2.5)*pow(I,ld-lp)*pow(-1.,l)*
							exp_delta_coulomb_f[lp]*exp_delta_coulomb_i[ld]*sqrt(2.*ld+1.))/(parm->k_Aa*k_p*sqrt(2.*l+1.));
					fl->energia=parm->energia_cm;
					fl->l=ld;
					fl->spin=0.;
					fl->j=ld;

					S[l]=GeneraDWspin(fl,&(parm->pot_opt[indx_ingreso]),parm->T_carga*parm->P_carga,parm->mu_Aa,
							parm->radio,parm->puntos,parm->matching_radio,&fp1);
					for(n=0;n<dim1->num_puntos;n++){
						rn=(dim1->a)+((dim1->b)-(dim1->a))*((dim1->puntos[n])+1.)/2.;
						for(m=0;m<=lp;m++){
							rhom[m]=0.;
						}
						rAn=km*rn;
						dim3->a=rAn-parm->r_A2max;
						dim3->b=rAn+parm->r_A2max;
						if(dim3->a<0.) dim3->a=0.;
						if(dim3->b>parm->radio) dim3->b=parm->radio-1.;
						GaussLegendre(dim3->puntos,dim3->pesos,dim3->num_puntos);
						SourcePrior2(rhom,nonm,fl,gl_up,gl_down,st,v_down,optico,core,l,rn,parm,dim3,dim2);
						for(m=0;m<=lp;m++){
							rho[n][l][m][lp]+=(redfac*rhofac*ClebsGordan(lp,-m,ld,0,l,-m)*rhom[0]);
							if(parm->prior==1) non[n][l][m][lp]+=(rhofac*ClebsGordan(lp,-m,ld,0,l,-m)*nonm[0]*rn);
						}
					}
				}
				dim1->a=parm->r_Ccmin;
				dim1->b=parm->r_Ccmax;
				if(energia_trans>0.) ElasticBreakup(Teb,rho,Ecm,v_up,v_down,dim1,parm,l,lp,k_n);
				for(n=0;n<dim1->num_puntos;n++){
					rn= (dim1->a)+((dim1->b)-(dim1->a))*((dim1->puntos[n])+1.)/2.;
					for(m=0;m<=lp;m++){
						phim[m]=0.;
					}
						NeutronWave(phim,rho,&(funcion_regular_up[0]),&(funcion_irregular_up[0]),dim1,parm,rn,l,lp,ld,wronskiano_up);
					for(m=0;m<=lp;m++){
					  	phi_up[n][l][m][lp]=((l+1.)/sqrt((l+1.)*(l+1.)+l*l))*phim[m];
					  	//phi_up[n][l][m][lp]=phim[m];
					}
					    NeutronWave(phim,rho,&(funcion_regular_down[1]),&(funcion_irregular_down[1]),dim1,parm,rn,l,lp,ld,wronskiano_down);
					for(m=0;m<=lp;m++){
					  phi_down[n][l][m][lp]=(l/sqrt((l+1.)*(l+1.)+l*l))*phim[m];
					  //phi_down[n][l][m][lp]=phim[m];
					}
				}
			}
			inc_break[l]=0.;
			elastic_break[l]=0.;
			inc_break_lmenos[l]=0.;
			inc_break_lmas[l]=rhoE*escala*sigma_const*AbsorcionPrior(direct,non_orth,cross_term,v_up,
					phi_up,non,dim1,l,parm->lmax,r_F);
			inc_break[l]=inc_break_lmas[l];
			inc_break_lmenos[l]=rhoE*escala*sigma_const*AbsorcionPrior(direct,non_orth,cross_term,v_down,
					phi_down,non,dim1,l,parm->lmax,r_F);

			inc_break[l]+=inc_break_lmenos[l];
			if(energia_trans>0.) elastic_break[l]=rhoE*rhoE_n*escala*sigma_const*PI*ElasticBreakupCross(Teb,l,parm->lmax);
			cross_total+=inc_break[l];
			cross_total_elasticb+=elastic_break[l];
			cross_total_breakup+=total_break[l];
			cout<<" NEB cross section: "<<inc_break[l]<<endl<<endl;
			cout<<" EB cross section: "<<elastic_break[l]<<endl<<endl;
			fp9<<"  "<<inc_break[l]<<"  "<<elastic_break[l]<<"  ";
			misc1<<l<<"  "<<inc_break[l]<<"  "<<elastic_break[l]<<endl;
		}
		TalysInput(inc_break_lmenos,inc_break_lmas,energia_trans,parm,&fp3,&fp4,&fp7,parm->J_A);
		cout<<"NEB cross section:  "<<cross_total<<"   EB cross section:  "<<cross_total_elasticb<<endl;
		fp9<<cross_total<<"  "<<cross_total_elasticb<<
				"  "<<cross_total+cross_total_elasticb<<endl;
		cross_total=0.;
		cross_total_elasticb=0.;
		for(l=0;l<parm->lmax;l++)
		{
			inc_break_lmenos[l]=0.;
			inc_break_lmas[l]=0.;
		}
  		cout<<"computing angular differential cross section"<<endl;
  		if(parm->capture_angular==1)
  		{
  			for(n=0;n<parm->cross_puntos;n++)
  			{
  				theta=PI*double(n)/double(parm->cross_puntos);
  				direct[0]=0.;
  				non_orth[0]=0.;
  				cross_term[0]=0.;
  				if((theta>=PI*parm->angle0/180.)&&(theta<=PI*parm->angle1/180.))
  				{
  					cross=AbsorcionAngular(v_up,phi_up,non,dim1,parm,theta,
  							direct,non_orth,cross_term,cross_up);
  					cross+=AbsorcionAngular(v_down,phi_down,non,dim1,parm,theta,
  							direct,non_orth,cross_term,cross_down);
  					elastic_cross=ElasticBreakupAngular(Teb,parm->lmax,theta);
  					cross_total+=sigma_const*escala*rhoE*cross*sin(theta)*2.*PI*PI/double(parm->cross_puntos);
  					cross_total_elasticb+=rhoE*rhoE_n*escala*sigma_const*PI*elastic_cross*sin(theta)*2.*PI*PI/double(parm->cross_puntos);
  					fp10<<theta*180./PI<<"  "<<sigma_const*escala*rhoE*cross<<
  							"  "<<rhoE*rhoE_n*escala*sigma_const*PI*elastic_cross<<"  "<<
							sigma_const*escala*rhoE*(cross)+(rhoE*rhoE_n*escala*sigma_const*PI*elastic_cross)<<endl;
  				}
  			}
  		}
		cout<<"NEB cross section:  "<<cross_total<<"   EB cross section  :  "<<cross_total_elasticb<<endl;
	}
	delete[] funcion_regular_up;
	delete[] funcion_irregular_up;
	delete[] funcion_regular_down;
	delete[] funcion_irregular_down;
	delete[] S;
	delete[] rho;
	delete[] rhom;
	delete[] non;
	delete[] nonm;
	delete[] phi_up;
	delete[] phi_down;
	delete[] phim;
	delete[] wf;
	delete[] rho_test;
	delete[] non_test;
	delete[] phi_test;
	delete[] phi_res;
	delete[] cross_up;
	delete[] cross_down;
	delete[] total_break;
	delete[] GreenFunction;
}
void AmplitudeCaptureCC(struct parametros* parm)
{

  parametros_integral *dim1=new parametros_integral;
  parametros_integral *dim2=new parametros_integral;
  parametros_integral *dim3=new parametros_integral;
  parametros_integral *dim4=new parametros_integral;
  complejo* exp_delta_coulomb_i=new complejo[parm->lmax];
  complejo* exp_delta_coulomb_f=new complejo[parm->lmax];
  double eta_f=parm->Z_a*parm->Z_A*E2HC*parm->mu_Bb*AMU/(HC*parm->k_Bb);
  double eta_i=parm->eta;
  double step,rn,energia_out,energia_trans,k_p,k_n,cross,elastic_cross,
    theta,costheta,D0,rhoE,sigma_const,escala,r_source,velocidad,
    cross_total,cross_total_elasticb,redfac,r_F,absorcion,e_res,rhoE_n,N_A,
    carga_out,carga_trans,km,rAn,Ecm,Ecm_out,cross_total_breakup,Ecmmax;
  distorted_wave* fl=new distorted_wave;
  distorted_wave* gl_up=new distorted_wave;
  distorted_wave* gl_down=new distorted_wave;
  potencial_optico *optico=new potencial_optico;
  potencial_optico *core=new potencial_optico;
  potencial_optico* v_up=new potencial_optico[1];
  potencial_optico* v_down=new potencial_optico[1];
  potencial_optico* vp_up=new potencial_optico[1];
  potencial_optico* vp_down=new potencial_optico[1];
  potencial_optico* pot_dumb=new potencial_optico;
  estado* st=new estado;
  estado* st_fin=new estado;
  ofstream fp1("dw_out1trans.txt");
  ofstream fp2("dw_in1trans.txt");
  ofstream fp3;
  fp3.open("talys1.txt");
  ofstream fp4;
  fp4.open("talys2.txt");
  ofstream fp5;
  fp5.open("talys_angular1.txt");
  ofstream fp6;
  fp6.open("talys_angular2.txt");
  ofstream fp7;
  fp7.open("SpinParity.txt");
  ofstream fp8;
  fp8.open("Jutta_angular.txt");
  ofstream fp9;
  fp9.open("dsdE.txt");
  ofstream fp10;
  fp10.open("dsdEdO.txt");
  complejo* S=new complejo[parm->lmax];
  complejo**** rho=tensor4_cmpx(parm->rCc_puntos,parm->lmax,parm->lmax+1,parm->lmax);
  complejo* rhom=new complejo[parm->lmax+1];
  complejo**** non=tensor4_cmpx(parm->rCc_puntos,parm->lmax,parm->lmax+1,parm->lmax);
  complejo**** dumb=tensor4_cmpx(parm->rCc_puntos,parm->lmax,parm->lmax+1,parm->lmax);
  complejo* nonm=new complejo[parm->lmax+1];
  complejo**** phi_up=tensor4_cmpx(parm->rCc_puntos,parm->lmax,parm->lmax+1,parm->lmax);
  complejo**** phi_down=tensor4_cmpx(parm->rCc_puntos,parm->lmax,parm->lmax+1,parm->lmax);
  complejo*** Teb=tensor_cmpx(parm->lmax,parm->lmax+1,parm->lmax);
  complejo** GreenFunction=matriz_cmpx(300,300);
  complejo** NLpot=matriz_cmpx(300,300);
  complejo* phi_res=new complejo [parm->puntos];
  complejo* phim=new complejo[parm->lmax+1];
  complejo* localgf=new complejo[1000];
  complejo* localpot=new complejo[1000];
  double* rg=new double[300];
  complejo pot_p;
  complejo pot_n;
  //pot opt_in(parm->dompot_n);
  //pot opt_out(parm->dompot_p);
  double* inc_break=new double[parm->lmax+1];
  double* inc_break_lmas=new double[parm->lmax+1];
  double* inc_break_lmenos=new double[parm->lmax+1];
  double* cross_up=new double[parm->lmax+1];
  double* cross_down=new double[parm->lmax+1];
  double* elastic_break=new double[parm->lmax+1];
  double* direct=new double[1];
  double* non_orth=new double[1];
  double* cross_term=new double[1];
  double** Al=matriz_dbl(2*parm->lmax+1,parm->lmax);
  int l,lp,dj,ld,indx_salida,indx_ingreso,indx_core,indx_neutron_target,indx_st,n,la,m,len,flag,n1,puntos_r,
    flagGF,flagpot;
  complejo rhofac,ampli,wronskiano,wronskiano_up,wronskiano_down,fl_int,gl_int,fl_source,gl_source,
    st_source,st_int,lorentz,is_pot_im,is_pot_im_out;
  dim1->a=parm->r_Ccmin;
  dim1->b=parm->r_Ccmax;
  dim1->num_puntos=parm->rCc_puntos;
  dim2->a=0.;
  dim2->b=PI;
  dim2->num_puntos=parm->theta_puntos;
  dim3->a=parm->r_A2min;
  dim3->b=parm->r_A2max;
  dim3->num_puntos=parm->rA2_puntos;
  dim4->a=parm->r_A2min;
  dim4->b=parm->r_A2max;
  dim4->num_puntos=parm->rA2_puntos;
  GaussLegendre(dim1->puntos,dim1->pesos,dim1->num_puntos);
  GaussLegendre(dim2->puntos,dim2->pesos,dim2->num_puntos);
  GaussLegendre(dim3->puntos,dim3->pesos,dim3->num_puntos);
  GaussLegendre(dim4->puntos,dim4->pesos,dim4->num_puntos);
  D0=10.;
  redfac=2.*AMU/(HC*HC);
  cout<<"Mass of projectile: "<<parm->P_masa<<endl;
  cout<<"Mass of target: "<<parm->T_masa<<endl;
  cout<<"Mass of residual nucleus: "<<parm->res_masa<<endl;
  carga_trans=parm->res_carga-parm->T_carga;
  cout<<"Charge of absorbed cluster: "<<carga_trans<<endl;
  carga_out=parm->P_carga-carga_trans;
  cout<<"Charge of emitted cluster: "<<carga_out<<endl;
  cout<<"Mass of absorbed cluster: "<<parm->n1_masa<<endl;
  cout<<"Mass of detected cluster: "<<parm->m_b<<endl;
  km=(parm->m_A+1.)/parm->m_A;
  if(parm->koning_delaroche==2) cout<<"******************************************************************"<<endl<<
      			              "***** Reading Coupled Cluster Green's function and potential *****"<<endl<<
				      "******************************************************************"<<endl;
  /*Selecciona los potenciales opticos en los distintos canales*/
  for (n=0;n<parm->num_opt;n++)
    {
      if(parm->optico_ingreso==parm->pot_opt[n].id) indx_ingreso=n;
      if(parm->optico_salida==parm->pot_opt[n].id) indx_salida=n;
      if(parm->optico_intermedio==parm->pot_opt[n].id) indx_neutron_target=n;
      if(parm->core_pot==parm->pot_opt[n].id) indx_core=n;
      //if(parm->pot_transfer==parm->pot_opt[n].id) v_up=&(parm->pot_opt[n]);
      //	if(parm->scatt_pot==parm->pot_opt[n].id) v_down=&(parm->pot_opt[n]);
    }
  //GeneraPotencialOptico(parm,v_up,1.,parm->m_A);
  //GeneraPotencialOptico(parm,v_down,1.,parm->m_A);
  v_up->puntos=parm->puntos;
  v_down->puntos=parm->puntos;
  vp_up->puntos=parm->puntos;
  vp_down->puntos=parm->puntos;
  v_up->radio=parm->radio;
  v_down->radio=parm->radio;
  vp_up->radio=parm->radio;
  vp_down->radio=parm->radio;
  v_up->radio_coul=5.;
  v_down->radio_coul=5.;
  vp_up->radio_coul=5.;
  vp_down->radio_coul=5.;
  v_up->Vso=5.;
  v_down->Vso=5.;
  vp_up->Vso=5.;
  vp_down->Vso=5.;
  v_up->radioso=5.;
  v_down->radioso=5.;
  vp_up->radioso=5.;
  vp_down->radioso=5.;
  v_up->aso=0.65;
  v_down->aso=0.65;
  vp_up->aso=0.65;
  vp_down->aso=0.65;
  cout<<"Center of mass energy: "<<parm->energia_cm<<endl;
  cout<<"initial momentum: "<<parm->k_Aa<<endl;
  cout<<"final momentum: "<<parm->k_Bb<<endl;
  cout<<"Initial reduced mass: "<<parm->mu_Aa<<endl;
  cout<<"final reduced mass: "<<parm->mu_Bb<<endl;
  cout<<"Absorbed cluster reduced mass: "<<parm->m_A/(parm->m_A+1.)<<endl;
  /*Calculo de las amplitudes de transferencia**************************************************************************/
  for(n=0;n<parm->num_st;n++)
    {
      if (parm->a_estados[0] == parm->st[n].id) st= &(parm->st[n]);
      if (parm->B_estados[0] == parm->st[n].id) st_fin= &(parm->st[n]);
    }
  step=double(parm->radio/parm->puntos);
  for(la=0;la<parm->lmax;la++)
    {
      exp_delta_coulomb_i[la]=exp(I*(deltac(la,eta_i)));
      exp_delta_coulomb_f[la]=exp(I*(deltac(la,eta_f)));
    }
  velocidad=C*sqrt(2*parm->energia_lab/(2.*AMU));
  sigma_const=2.*parm->mu_Aa*AMU/(HC*HC*parm->k_Aa);
  len=strlen(parm->unidades);
  if(!strncmp(parm->unidades,"milib",len)) flag=1;
  if(!strncmp(parm->unidades,"fm2",len)) flag=2;
  if(!strncmp(parm->unidades,"b",len)) flag=3;
  if(!strncmp(parm->unidades,"microb",len)) flag=4;
  switch(flag)
    {
    case 1:
      escala=10.;
      cout<<"Cross section in milibarn"<<endl;
      break;
    case 2:
      escala=1.;
      cout<<"Cross section in fm^2"<<endl;
      break;
    case 3:
      escala=0.01;
      cout<<"Cross section in barn"<<endl;
      break;
    case 4:
      escala=10000.;
      cout<<"Cross section in microbarn"<<endl;
      break;
    default:
      Error("Unkown units for cros section");
      break;
    }
  r_F=1000.;
  e_res=st_fin->energia;
  puntos_r=200;
  energia_out=parm->enerange_max;
  Ecm_out=(parm->T_masa)*energia_out/(parm->n1_masa+(parm->T_masa));
  Ecmmax=parm->energia_cm-(((parm->T_masa)*parm->enerange_min/(parm->n1_masa+(parm->T_masa))))-2.2245;
  Ecm=parm->energia_cm-Ecm_out-2.2245;
  for(;;)
    {
      l=1;
      dj=3;      
      flagGF=ReadGF("/home/gregory/DOM/localdom/potential_Gregory/BIG_gf_reg_mesh",GreenFunction,rg,
		    puntos_r,&Ecm,l,Ecmmax,dj);
      flagpot=ReadNLpot("/home/gregory/DOM/localdom/potential_Gregory/BIG_optical_pot",NLpot,rg,
		    puntos_r,Ecm,l,dj);
      //Localize(NLpot,rg,puntos_r,localpot,dim1);
      misc2<<Ecm<<"  "<<real(GreenFunction[20][20])<<"  "<<imag(GreenFunction[20][20])<<endl;
      //Localize(GreenFunction,rg,puntos_r,localpot,dim1);
      //exit(0);
      if(flagGF==0 || flagpot==0)
	{
	  cout<<"Exiting loop on Ecm="<<Ecm<<endl;
	  break;
	}
      cout<<"flagGF: "<<flagGF<<"   flagpot: "<<flagpot<<endl;
      Ecm_out=parm->energia_cm-Ecm-2.2245;
      energia_out=(parm->n1_masa+(parm->T_masa))*Ecm_out/(parm->T_masa);
      energia_trans=(parm->n1_masa+parm->T_masa)*Ecm/(parm->T_masa);
      cout<<"Energy of detected cluster: "<<energia_out<<endl<<"Energy of absorbed cluster: "
	  <<energia_trans<<endl<<"CM energy of absorbed cluster-target system: "
	  <<Ecm<<endl<<"CM energy of detected cluster-target system: "
	  <<Ecm_out<<endl;
      cout<<"Conditions: "<<(Ecm<Ecmmax)<<"  "<<(flagGF)<<endl;
      fp9<<energia_out<<"  "<<energia_trans<<"  ";
      misc1<<"& Energy of detected cluster: "<<energia_out<<"    Energy of absorbed cluster: "<<energia_trans<<endl;
      //		misc2<<endl<<"*********************  Ep= "<<energia_out<<" ****************************"<<endl;
      k_n=sqrt(2.*parm->n1_masa*AMU*Ecm)/HC;
      k_p=sqrt(2.*parm->m_b*AMU*Ecm_out)/HC;
      rhoE=parm->m_b*AMU*k_p/(8.*PI*PI*PI*HC*HC);
      rhoE_n=parm->n1_masa*AMU*k_n/(8.*PI*PI*PI*HC*HC);
      eta_f=carga_out*parm->res_carga*E2HC*(parm->m_b*parm->T_masa/(parm->m_b+parm->T_masa))*AMU/(HC*k_p);
      cross_total=0.;
      cross_total_elasticb=0.;
      //exit(0);
      for(l=parm->lmin;l<parm->ltransfer;l++)
	{
	  cout<<"L: "<<l<<endl;
	  for(n=0;n<parm->puntos;n++)
	    {
	      rn=step*(n+1.);
	      KoningDelaroche(energia_out,parm->T_N,parm->T_carga,rn,&pot_p,
			      &pot_n,l,l+0.5,vp_down,pot_dumb);
	      vp_up->r[n]=rn;
	      vp_up->pot[n]=pot_p;
	      if(l>0) KoningDelaroche(energia_out,parm->T_N,parm->T_carga,rn,&pot_p,
				      &pot_n,l,l-0.5,vp_down,pot_dumb);
	      if(l==0) KoningDelaroche(energia_out,parm->T_N,parm->T_carga,rn,&pot_p,
				       &pot_n,l,0.5,vp_down,pot_dumb);
	      vp_down->r[n]=rn;
	      vp_down->pot[n]=pot_p;
	      KoningDelaroche(energia_trans,parm->T_N,parm->T_carga,rn,&pot_p,
			      &pot_n,l,l+0.5,pot_dumb,v_down);
	      v_up->r[n]=rn;
	      v_up->pot[n]=pot_n;
	      v_down->pot[n]=pot_n;
	      if(l>0) KoningDelaroche(energia_trans,parm->T_N,parm->T_carga,rn,&pot_p,
				      &pot_n,l,l-0.5,pot_dumb,v_down);
	      if(l==0) KoningDelaroche(energia_out,parm->T_N,parm->T_carga,rn,&pot_p,
				       &pot_n,l,0.5,pot_dumb,v_down);
	      v_down->r[n]=rn;
	      v_down->pot[n]=pot_n;
	    }
	  for(lp=0;lp<parm->lmax;lp++)
	    {
	      if(parm->remnant==1 && parm->prior==1) {
		GeneraRemnant(optico,core,&parm->pot_opt[indx_ingreso],vp_down,parm->T_carga*parm->P_carga,
			      0.,0,0,parm->mu_Aa,parm->m_b);
	      }
	      gl_up->energia=Ecm_out;
	      gl_up->l=lp;
	      gl_up->spin=parm->n_spin;
	      gl_up->j=lp+parm->n_spin;
	      GeneraDWspin(gl_up,vp_up,0.,parm->m_b*parm->res_masa/(parm->m_b+parm->res_masa),
			   parm->radio,parm->puntos,parm->matching_radio,&fp2);
	      gl_down->energia=Ecm_out;
	      gl_down->l=lp;
	      gl_down->spin=parm->n_spin;
	      gl_down->j=lp-parm->n_spin;
	      if(lp==0) gl_down->j=lp;
	      GeneraDWspin(gl_down,vp_down,0.,parm->m_b*parm->res_masa/(parm->m_b+parm->res_masa),
			   parm->radio,parm->puntos,parm->matching_radio,&fp2);
	      for(n=0;n<dim1->num_puntos;n++){
		for(m=0;m<=lp;m++){
		  rho[n][l][m][lp]=0.;
		  non[n][l][m][lp]=0.;
		}
	      }
	      exp_delta_coulomb_f[lp]=exp(I*(deltac(lp,eta_f)));
	      for(ld=abs(l-lp);(ld<=l+lp)&&(ld<parm->lmax);ld++)
		{
		  rhofac=(16.*pow(PI,2.5)*pow(I,ld-lp)*pow(-1.,l)*
			  exp_delta_coulomb_f[lp]*exp_delta_coulomb_i[ld]*sqrt(2.*ld+1.))/(parm->k_Aa*k_p*sqrt(2.*l+1.));
		  fl->energia=parm->energia_cm;
		  fl->l=ld;
		  fl->spin=0.;
		  fl->j=ld;

		  S[l]=GeneraDWspin(fl,&(parm->pot_opt[indx_ingreso]),parm->T_carga*parm->P_carga,parm->mu_Aa,
				    parm->radio,parm->puntos,parm->matching_radio,&fp1);
		  for(n=0;n<dim1->num_puntos;n++){
		    rn=(dim1->a)+((dim1->b)-(dim1->a))*((dim1->puntos[n])+1.)/2.;
		    for(m=0;m<=lp;m++){
		      rhom[m]=0.;
		    }
		    rAn=km*rn;
		    dim3->a=rAn-parm->r_A2max;
		    dim3->b=rAn+parm->r_A2max;
		    if(dim3->a<0.) dim3->a=0.;
		    if(dim3->b>parm->radio) dim3->b=parm->radio-1.;
		    GaussLegendre(dim3->puntos,dim3->pesos,dim3->num_puntos);
		    //SourcePrior2(rhom,nonm,fl,gl_up,gl_down,st,v_down,optico,core,l,rn,parm,dim3,dim2);
		    SourceNL(rhom,nonm,fl,gl_up,gl_down,st,NLpot,rg,puntos_r,optico,core,l,rn,parm,dim3,dim2);
		    for(m=0;m<=lp;m++){
		      rho[n][l][m][lp]+=(redfac*rhofac*ClebsGordan(lp,-m,ld,0,l,-m)*rhom[0]);
		      if(parm->prior==1) non[n][l][m][lp]+=(rhofac*ClebsGordan(lp,-m,ld,0,l,-m)*nonm[0]*rn);
		    }
		    //						misc2<<rn<<"  "<<real(rho[n][l][0][lp])<<"  "<<imag(rho[n][l][0][lp])<<endl;
		  }
		  //					exit(0);
		}
	      dim1->a=parm->r_Ccmin;
	      dim1->b=parm->r_Ccmax;
	      if(energia_trans>0.) ElasticBreakup(Teb,rho,Ecm,v_up,v_down,dim1,parm,l,lp,k_n);
	      for(n=0;n<dim1->num_puntos;n++){
		rn= (dim1->a)+((dim1->b)-(dim1->a))*((dim1->puntos[n])+1.)/2.;
		for(m=0;m<=lp;m++){
		  phim[m]=0.;
		}
		NeutronWaveGF(phim,rho,GreenFunction,rg,puntos_r,dim1,parm,rn,l,lp,ld,k_n);
		for(m=0;m<=lp;m++){
		  phi_up[n][l][m][lp]=phim[m];
		}
		if ((Ecm>-4.1)&&(Ecm<-3.9))  misc4<<rn<<"  "<<real(phi_up[n][l][0][lp])
						  <<"  "<<imag(phi_up[n][l][0][lp])<<"  "<<abs(phi_up[n][l][0][lp])<<endl;
	      }
	      //exit(0);
	    }
	  inc_break[l]=0.;
	  elastic_break[l]=0.;
	  inc_break_lmenos[l]=0.;
	  // inc_break_lmas[l]=rhoE*escala*sigma_const*AbsorcionPrior(direct,non_orth,cross_term,v_up,
	  // 							   phi_up,non,dim1,l,parm->lmax,r_F);
	  inc_break_lmas[l]=rhoE*escala*sigma_const*AbsorcionNL(NLpot,phi_up,non,dim1,l,parm->lmax,rg,puntos_r);
	  inc_break[l]=inc_break_lmas[l];
	  //			inc_break_lmenos[l]=rhoE*escala*sigma_const*AbsorcionPrior(direct,non_orth,cross_term,v_down,
	  //					phi_down,dumb,dim1,l,parm->lmax,r_F);

	  //			inc_break[l]+=inc_break_lmenos[l];
	  if(energia_trans>0.) elastic_break[l]=rhoE*rhoE_n*escala*sigma_const*PI*ElasticBreakupCross(Teb,l,parm->lmax);
	  cross_total+=inc_break[l];
	  cross_total_elasticb+=elastic_break[l];
	  cout<<" NEB cross section: "<<inc_break[l]<<endl<<endl;
	  cout<<" EB cross section: "<<elastic_break[l]<<endl<<endl;
	  fp9<<"  "<<inc_break[l]<<"  "<<elastic_break[l]<<"  ";
	  misc1<<l<<"  "<<inc_break[l]<<"  "<<elastic_break[l]<<endl;
	}
      TalysInput(inc_break_lmenos,inc_break_lmas,energia_trans,parm,&fp3,&fp4,&fp7,parm->J_A);
      cout<<"NEB cross section:  "<<cross_total<<"   EB cross section:  "<<cross_total_elasticb<<endl;
      fp9<<cross_total<<"  "<<cross_total_elasticb<<
	"  "<<cross_total+cross_total_elasticb<<endl;
      cross_total=0.;
      cross_total_elasticb=0.;
      for(l=0;l<parm->lmax;l++)
	{
	  inc_break_lmenos[l]=0.;
	  inc_break_lmas[l]=0.;
	}
      if(parm->capture_angular==1)
	{
	  cout<<"computing angular differential cross section"<<endl;
	  for(n=0;n<parm->cross_puntos;n++)
	    {
	      theta=PI*double(n)/double(parm->cross_puntos);
	      direct[0]=0.;
	      non_orth[0]=0.;
	      cross_term[0]=0.;
	      cross=0.;
	      if((theta>=PI*parm->angle0/180.)&&(theta<=PI*parm->angle1/180.))
		{
		  cross=AbsorcionAngular(v_up,phi_up,non,dim1,parm,theta,
					 direct,non_orth,cross_term,cross_up);
		  //			cross+=AbsorcionAngular(v_down,phi_down,non,dim1,parm,theta,
		  //					direct,non_orth,cross_term,cross_down);
		  elastic_cross=ElasticBreakupAngular(Teb,parm->lmax,theta);
		  cross_total+=sigma_const*escala*rhoE*cross*sin(theta)*2.*PI*PI/double(parm->cross_puntos);
		  cross_total_elasticb+=rhoE*rhoE_n*escala*sigma_const*PI*elastic_cross*sin(theta)*2.*PI*PI/double(parm->cross_puntos);
		  fp10<<theta*180./PI<<"  "<<sigma_const*escala*rhoE*cross<<
		    "  "<<rhoE*rhoE_n*escala*sigma_const*PI*elastic_cross<<"  "<<
		    sigma_const*escala*rhoE*(cross)+(rhoE*rhoE_n*escala*sigma_const*PI*elastic_cross)<<endl;
		}
	    }
	  cout<<"NEB cross section:  "<<cross_total<<"   EB cross section:  "<<cross_total_elasticb<<endl;
	}
      cout<<"Conditions: "<<(Ecm<Ecmmax)<<"  "<<(flagGF)<<endl;
    }
  cout<<"Out of loop"<<endl;
  delete[] S;
  delete[] rho;
  delete[] rhom;
  delete[] non;
  delete[] nonm;
  delete[] phi_up;
  delete[] phi_down;
  delete[] phim;
  delete[] localpot;
  delete[] localgf;
  delete[] phi_res;
  delete[] cross_up;
  delete[] cross_down;
  delete[] GreenFunction;
}
void Source(complejo* rho,distorted_wave* f,distorted_wave* g,estado* u,potencial* v,int l,
		double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2)
{
	int n1,n2,n3,m,lp,ld;
	double rAp,theta,rd,rdx,rdz,rpnx,rpnz,rpn,vpn,seno,coseno,coseno_d,
	          theta_Bn,seno_Bn,coseno_Bn,km;
	complejo fl,gl,ud;
	complejo* suma=new complejo[l+1];
	ld=f->l;
	lp=g->l;
	km=parm->m_A/(parm->m_A+1.);
	for(m=0;m<=l;m++){
		suma[m]=0.;
	}
	for (n1 =0;n1<dim1->num_puntos; n1++) {
		rAp = (dim1->a)+((dim1->b)-(dim1->a))*((dim1->puntos[n1])+1.)/2.;
		fl=interpola_cmpx(f->wf,f->r,rAp,f->puntos);
		for (n2=0;n2<dim2->num_puntos;n2++) {
			theta=(dim2->a)+((dim2->b)-(dim2->a))*((dim2->puntos[n2])+1.)/2.;
			seno=sin(theta);
			coseno=cos(theta);
			for (n3=0;n3<dim2->num_puntos;n3++) {
				theta_Bn=(dim2->a)+((dim2->b)-(dim2->a))*((dim2->puntos[n3])+1.)/2.;
				seno_Bn=sin(theta_Bn);
				coseno_Bn=cos(theta_Bn);
				rdx=0.5*(rAp*seno+km*rBn*seno_Bn);
				rdz=0.5*(rAp*coseno+km*rBn*coseno_Bn);
				rd=sqrt(rdx*rdx+rdz*rdz);
				rpnx=rBn*km*seno_Bn-rAp*seno;
				rpnz=rBn*km*coseno_Bn-rAp*coseno;
				rpn=sqrt(rpnx*rpnx+rpnz*rpnz);
				coseno_d=rdz/rd;
				gl=interpola_cmpx(g->wf,g->r,rd,g->puntos);
				ud=interpola_cmpx(u->wf,u->r,rpn,u->puntos);
				vpn=interpola_dbl(v->pot,v->r,rpn,v->puntos);
				for(m=0;m<=l;m++){
					if(m<=lp) suma[m]+=rAp*seno*seno_Bn*fl*gl*ud*vpn*FuncionAngular(lp,ld,l,m,l,m,coseno,coseno_d,coseno_Bn)*
							dim1->pesos[n1]*dim2->pesos[n2]*dim2->pesos[n3]/(rd);
				}
			}
		}
	}
	for(m=0;m<=l;m++){
		rho[m]=suma[m]*((dim1->b)-(dim1->a))*((dim2->b)-(dim2->a))*((dim2->b)-(dim2->a))/8.;
	}
	delete[] suma;
}
void Source2(complejo* rho,distorted_wave* f,distorted_wave* g,estado* u,potencial* v,int l,
		double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2)
{
	int n1,n2,n3,m,lp,ld;
	double rAp,theta,rd,rdx,rdz,rpnx,rpnz,rpn,vpn,seno,coseno,coseno_d,
	          theta_Bn,seno_Bn,coseno_Bn,km;
	complejo fl,gl,ud;
	complejo* suma=new complejo[l+1];
	ld=f->l;
	lp=g->l;
	km=parm->m_A/(parm->m_A+1.);
	for(m=0;m<=l;m++){
		suma[m]=0.;
	}
	for (n1 =0;n1<dim1->num_puntos; n1++) {
		rAp = (dim1->a)+((dim1->b)-(dim1->a))*((dim1->puntos[n1])+1.)/2.;
		fl=interpola_cmpx(f->wf,f->r,rAp,f->puntos);
		for (n2=0;n2<dim2->num_puntos;n2++) {
			theta=(dim2->a)+((dim2->b)-(dim2->a))*((dim2->puntos[n2])+1.)/2.;
			seno=sin(theta);
			coseno=cos(theta);
			for (n3=0;n3<dim2->num_puntos;n3++) {
				theta_Bn=(dim2->a)+((dim2->b)-(dim2->a))*((dim2->puntos[n3])+1.)/2.;
				seno_Bn=sin(theta_Bn);
				coseno_Bn=cos(theta_Bn);
				rdx=0.5*(rAp*seno+km*rBn*seno_Bn);
				rdz=0.5*(rAp*coseno+km*rBn*coseno_Bn);
				rd=sqrt(rdx*rdx+rdz*rdz);
				rpnx=rBn*km*seno_Bn-rAp*seno;
				rpnz=rBn*km*coseno_Bn-rAp*coseno;
				rpn=sqrt(rpnx*rpnx+rpnz*rpnz);
				coseno_d=rdz/rd;
				gl=interpola_cmpx(g->wf,g->r,rd,g->puntos);
				ud=interpola_cmpx(u->wf,u->r,rpn,u->puntos);
				vpn=interpola_dbl(v->pot,v->r,rpn,v->puntos);
				for(m=0;m<=l;m++){
					if(m<=lp) suma[m]+=rAp*seno*seno_Bn*fl*gl*ud*vpn*FuncionAngular(lp,ld,l,m,l,m,coseno,coseno_d,coseno_Bn)*
							dim1->pesos[n1]*dim2->pesos[n2]*dim2->pesos[n3]/(rd);
				}
			}
		}
	}
	for(m=0;m<=l;m++){
		rho[m]=suma[m]*((dim1->b)-(dim1->a))*((dim2->b)-(dim2->a))*((dim2->b)-(dim2->a))/8.;
	}
	delete[] suma;
}
void SourcePrior(complejo* rho,complejo* non,distorted_wave* f,distorted_wave* g,estado* u,potencial* v,int l,
		double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2)
{
	int n1,n2,n3,m,lp,ld;
	double rAp,theta,rd,rdx,rdz,rpnx,rpnz,rpn,vpn,seno,coseno,coseno_d,
	          theta_Bn,seno_Bn,coseno_Bn,km,rAn;
	complejo fl,gl,ud,coupling;
	complejo* suma=new complejo[l+1];
	complejo* sumanon=new complejo[l+1];
	ld=f->l;
	lp=g->l;
	km=(parm->m_A+1.)/parm->m_A;
	for(m=0;m<=l;m++){
		suma[m]=0.;
		sumanon[m]=0.;
	}
	rAn=km*rBn;
	vpn=interpola_dbl(v->pot,v->r,rAn,v->puntos);
	for (n1 =0;n1<dim1->num_puntos; n1++) {
		rAp = (dim1->a)+((dim1->b)-(dim1->a))*((dim1->puntos[n1])+1.)/2.;
		fl=interpola_cmpx(f->wf,f->r,rAp,f->puntos);
		for (n2=0;n2<dim2->num_puntos;n2++) {
			theta=(dim2->a)+((dim2->b)-(dim2->a))*((dim2->puntos[n2])+1.)/2.;
			seno=sin(theta);
			coseno=cos(theta);
			for (n3=0;n3<dim2->num_puntos;n3++) {
				theta_Bn=(dim2->a)+((dim2->b)-(dim2->a))*((dim2->puntos[n3])+1.)/2.;
				seno_Bn=sin(theta_Bn);
				coseno_Bn=cos(theta_Bn);
				rdx=0.5*(rAp*seno+km*rBn*seno_Bn);
				rdz=0.5*(rAp*coseno+km*rBn*coseno_Bn);
				rd=sqrt(rdx*rdx+rdz*rdz);
				rpnx=rBn*km*seno_Bn-rAp*seno;
				rpnz=rBn*km*coseno_Bn-rAp*coseno;
				rpn=sqrt(rpnx*rpnx+rpnz*rpnz);
				coseno_d=rdz/rd;
				gl=interpola_cmpx(g->wf,g->r,rd,g->puntos);
				ud=interpola_cmpx(u->wf,u->r,rpn,u->puntos);
				for(m=0;m<=l;m++){
					if(m<=lp) {
						coupling=FuncionAngular(lp,ld,l,m,l,m,coseno,coseno_d,coseno_Bn);
						suma[m]+=rAp*seno*seno_Bn*fl*gl*ud*vpn*coupling*
							dim1->pesos[n1]*dim2->pesos[n2]*dim2->pesos[n3]/(rd);
						sumanon[m]+=rAp*seno*seno_Bn*fl*gl*ud*coupling*
							dim1->pesos[n1]*dim2->pesos[n2]*dim2->pesos[n3]/(rd);
					}
				}
			}
		}
	}
	for(m=0;m<=l;m++){
		rho[m]=suma[m]*((dim1->b)-(dim1->a))*((dim2->b)-(dim2->a))*((dim2->b)-(dim2->a))/8.;
		non[m]=sumanon[m]*((dim1->b)-(dim1->a))*((dim2->b)-(dim2->a))*((dim2->b)-(dim2->a))/8.;
	}
//	misc1<<rBn<<"  "<<abs(rho[0])<<endl;
	delete[] suma;
	delete[] sumanon;
}
void SourcePrior2(complejo* rho,complejo* non,distorted_wave* f,distorted_wave* g_up,
		  distorted_wave* g_down,estado* u,potencial_optico* v,potencial_optico* optico,
	 potencial_optico* core,int l,double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2)
{
  int n1,n2,m,lp,ld;
  double rAp,theta,rd,rdx,rdz,rpnx,rpnz,rpn,seno,coseno,coseno_d,
    km,rAn,interruptor,rBp,rBpx,rBpz,upfrac,downfrac;
  complejo fl,gl_up,gl_down,ud,coupling,corepot,inpot,remnant,vpn,suma,sumanon,sumanonZR,N0;
  ld=f->l;
  lp=g_up->l;
  km=(parm->m_A+1.)/parm->m_A;
  suma=0.;
  sumanon=0.;
  sumanonZR=0.;
  rAn=km*rBn;
  N0=0.;
  upfrac=(lp+1.)/(sqrt((lp+1.)*(lp+1.)+lp*lp));
  downfrac=(lp/(sqrt((lp+1.)*(lp+1.)+lp*lp)));
  upfrac=1.;
  downfrac=0.;
  vpn=interpola_cmpx(v->pot,v->r,rAn,v->puntos);
  remnant=0.;
  for (n1 =0;n1<dim1->num_puntos; n1++) {
    rAp = (dim1->a)+((dim1->b)-(dim1->a))*((dim1->puntos[n1])+1.)/2.;
    if(parm->remnant==1) corepot=interpola_cmpx(core->pot,core->r,rAp,core->puntos);
    for (n2=0;n2<dim2->num_puntos;n2++) {
      theta=(dim2->a)+((dim2->b)-(dim2->a))*((dim2->puntos[n2])+1.)/2.;
      //			theta=PI/2.;
      seno=sin(theta);
      coseno=cos(theta);
      rdx=0.5*(rAp*seno);
      rdz=0.5*(rAp*coseno+rAn);
      rd=sqrt(rdx*rdx+rdz*rdz);
      if(parm->remnant==1) inpot=interpola_cmpx(optico->pot,optico->r,rd,optico->puntos);
      rpnx=-rAp*seno;
      rpnz=rAn-rAp*coseno;
      rpn=sqrt(rpnx*rpnx+rpnz*rpnz);
      coseno_d=rdz/rd;
      rBpx=rAp*seno;
      rBpz=(-1./parm->m_A)*rBn+rAp*coseno;
      rBp=sqrt(rBpx*rBpx+rBpz*rBpz);
      gl_up=interpola_cmpx(g_up->wf,g_up->r,rBp,g_up->puntos);
      gl_down=interpola_cmpx(g_down->wf,g_down->r,rBp,g_down->puntos);
      fl=interpola_cmpx(f->wf,f->r,rd,f->puntos);
      ud=interpola_cmpx(u->wf,u->r,rpn,u->puntos);
      coupling=FuncionAngular2(lp,ld,l,coseno,coseno_d);
      if (parm->remnant==1) remnant=inpot-corepot;
      suma+=rBp*seno*fl*(upfrac*gl_up+downfrac*gl_down)*ud*(vpn-remnant)*coupling*
	dim1->pesos[n1]*dim2->pesos[n2]/(rd);
      sumanon+=rBp*seno*fl*(upfrac*gl_up+downfrac*gl_down)*ud*coupling*
	dim1->pesos[n1]*dim2->pesos[n2]/(rd);
    }
  }
  rho[0]=suma*((dim1->b)-(dim1->a))*((dim2->b)-(dim2->a))/4.;
  non[0]=sumanon*((dim1->b)-(dim1->a))*((dim2->b)-(dim2->a))/4.;
}
void SourceNL(complejo* rho,complejo* non,distorted_wave* f,distorted_wave* g_up,
	      distorted_wave* g_down,estado* u,complejo** v,double* r,int puntos_r,potencial_optico* optico,
	 potencial_optico* core,int l,double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2)
{
  int n1,n2,n3,m,lp,ld;
  double rAp,theta,rd,rdx,rdz,rpnx,rpnz,rpn,seno,coseno,coseno_d,
    km,rAn,interruptor,rBp,rBpx,rBpz,upfrac,downfrac,rAnnl,rdnl,rdxnl,rdznl, coseno_dnl
    ,rpnnl,rpxnl,rpnznl;
  complejo fl,gl_up,gl_down,ud,coupling,corepot,inpot,remnant,vpn,suma,sumanon,sumanonZR,N0;
  ld=f->l;
  lp=g_up->l;
  km=(parm->m_A+1.)/parm->m_A;
  suma=0.;
  sumanon=0.;
  sumanonZR=0.;
  rAn=km*rBn;
  N0=0.;
  upfrac=(lp+1.)/(sqrt((lp+1.)*(lp+1.)+lp*lp));
  downfrac=(lp/(sqrt((lp+1.)*(lp+1.)+lp*lp)));
  upfrac=1.;
  downfrac=0.;
  remnant=0.;
  for (n1 =0;n1<dim1->num_puntos; n1++) {
    rAp = (dim1->a)+((dim1->b)-(dim1->a))*((dim1->puntos[n1])+1.)/2.;
    if(parm->remnant==1) corepot=interpola_cmpx(core->pot,core->r,rAp,core->puntos);
    for (n2=0;n2<dim2->num_puntos;n2++) {
      theta=(dim2->a)+((dim2->b)-(dim2->a))*((dim2->puntos[n2])+1.)/2.;
      seno=sin(theta);
      coseno=cos(theta);
      rdx=0.5*(rAp*seno);
      rdz=0.5*(rAp*coseno+rAn);
      rd=sqrt(rdx*rdx+rdz*rdz);
      if(parm->remnant==1) inpot=interpola_cmpx(optico->pot,optico->r,rd,optico->puntos);
      rpnx=-rAp*seno;
      rBpx=rAp*seno;
      rBpz=(-1./parm->m_A)*rBn+rAp*coseno;
      rBp=sqrt(rBpx*rBpx+rBpz*rBpz);
      gl_up=interpola_cmpx(g_up->wf,g_up->r,rBp,g_up->puntos);
      gl_down=interpola_cmpx(g_down->wf,g_down->r,rBp,g_down->puntos);
      for (n3=0;n3<dim1->num_puntos; n3++) {
	rAnnl = (dim1->a)+((dim1->b)-(dim1->a))*((dim1->puntos[n3])+1.)/2.;
	vpn=interpola2D_cmpx(v,r,r,rAn,rAnnl,puntos_r,puntos_r);
	rdznl=0.5*(rAp*coseno+rAnnl);
	rdnl=sqrt(rdx*rdx+rdznl*rdznl);
	fl=interpola_cmpx(f->wf,f->r,rdnl,f->puntos);
	rpnznl=rAnnl-rAp*coseno;
	rpnnl=sqrt(rpnx*rpnx+rpnznl*rpnznl); 
	ud=interpola_cmpx(u->wf,u->r,rpnnl,u->puntos);
	coseno_dnl=rdznl/rdnl;
	coupling=FuncionAngular2(lp,ld,l,coseno,coseno_d);
	if (parm->remnant==1) remnant=inpot-corepot;
	suma+=rBp*seno*fl*(upfrac*gl_up+downfrac*gl_down)*ud*(vpn-remnant)*coupling*
	  dim1->pesos[n1]*dim1->pesos[n3]*dim2->pesos[n2]/(rdnl);
	sumanon+=rBp*seno*fl*(upfrac*gl_up+downfrac*gl_down)*ud*coupling*
	  dim1->pesos[n1]*dim1->pesos[n3]*dim2->pesos[n2]/(rdnl);
      }
    }
  }
  rho[0]=suma*((dim1->b)-(dim1->a))*((dim1->b)-(dim1->a))*((dim2->b)-(dim2->a))/8.;
  non[0]=sumanon*((dim1->b)-(dim1->a))*((dim1->b)-(dim1->a))*((dim2->b)-(dim2->a))/8.;
}
void Localize(complejo** nlf,double* r,int puntos_r,complejo* lf,parametros_integral* dim)
{
  int n1,n2;
  double R,Rp;
  complejo nlfint;
  for (n1 =0;n1<puntos_r; n1++) {
    lf[n1]=0.;
    R=r[n1];
    for (n2 =0;n2<dim->num_puntos; n2++) {
      Rp= (dim->a)+((dim->b)-(dim->a))*((dim->puntos[n2])+1.)/2.;
      nlfint=interpola2D_cmpx(nlf,r,r,R,Rp,puntos_r,puntos_r);
      lf[n1]+=Rp*Rp*nlfint*dim->pesos[n2]*((dim->b)-(dim->a))/2.;
    }
    misc1<<R<<"  "<<real(lf[n1])<<"  "<<imag(lf[n1])<<endl;
    misc2<<R<<"  "<<real(nlf[n1][10])<<"  "<<imag(nlf[n1][10])<<endl;
  }
}

void TestIntegral(distorted_wave* f,distorted_wave* g,estado* u,potencial* v,int l,int m,int K,int M,
		double rBn,parametros* parm, parametros_integral* dim1,parametros_integral* dim2)
{
	int n1,n2,n3,n4,n5,lp,ld;
	double rAp,theta,rd,rdx,rdz,rdy,rpnx,rpnz,rpny,rpn,vpn,seno,coseno,coseno_d,
	          km,rAn,thetaBn,phi,phiBn,senoBn,cosenoBn,phi_d;
	complejo fl,gl,ud,coupling,I1,I2;
	ld=f->l;
	lp=g->l;
	km=(parm->m_A+1.)/parm->m_A;
	rAn=km*rBn;
	vpn=interpola_dbl(v->pot,v->r,rAn,v->puntos);
	I1=0.;
	rAp=3.;
	fl=interpola_cmpx(f->wf,f->r,rAp,f->puntos);
	for (n2=0;n2<dim2->num_puntos;n2++) {
		theta=(dim2->a)+((dim2->b)-(dim2->a))*((dim2->puntos[n2])+1.)/2.;
		seno=sin(theta);
		coseno=cos(theta);
		for (n3=0;n3<dim2->num_puntos;n3++) {
			thetaBn=(dim2->a)+((dim2->b)-(dim2->a))*((dim2->puntos[n3])+1.)/2.;
			senoBn=sin(thetaBn);
			cosenoBn=cos(thetaBn);
			for (n4=0;n4<dim2->num_puntos;n4++) {
				phi=2.*PI*((dim2->puntos[n4])+1.)/2.;
				for (n5=0;n5<dim2->num_puntos;n5++) {
					phiBn=2.*PI*((dim2->puntos[n5])+1.)/2.;
					rdx=0.5*(rAp*seno*cos(phi)+km*rBn*senoBn*cos(phiBn));
					rdy=0.5*(rAp*seno*sin(phi)+km*rBn*senoBn*sin(phiBn));
					rdz=0.5*(rAp*coseno+km*rBn*cosenoBn);
					rd=sqrt(rdx*rdx+rdy*rdy+rdz*rdz);
					rpnx=rBn*km*senoBn*cos(phiBn)-rAp*seno*cos(phi);
					rpny=rBn*km*senoBn*sin(phiBn)-rAp*seno*sin(phi);
					rpnz=rBn*km*cosenoBn-rAp*coseno;
					rpn=sqrt(rpnx*rpnx+rpny*rpny+rpnz*rpnz);
					coseno_d=rdz/rd;
					phi_d=atan2(rdy,rdx);
					gl=interpola_cmpx(g->wf,g->r,rd,g->puntos);
					ud=interpola_cmpx(u->wf,u->r,rpn,u->puntos);
					coupling=FuncionAngular3(lp,ld,K,M,coseno,coseno_d,phi,phi_d)*SphericalHarmonic(l,-m,cosenoBn,phiBn);
//					coupling=FuncionAngular3(1,1,2,0,coseno,coseno_d,phi,phi_d)*SphericalHarmonic(2,0,cosenoBn,phiBn);
					I1+=seno*senoBn*fl*gl*ud*vpn*coupling*
							dim2->pesos[n2]*dim2->pesos[n3]*dim2->pesos[n4]*dim2->pesos[n5]/(rd);
//					I1+=seno*senoBn*coupling*gl*ud*fl*dim2->pesos[n2]*dim2->pesos[n3]*dim2->pesos[n4]*dim2->pesos[n5];
				}
			}
		}
	}
	I1=I1*PI*PI*PI*PI/4.;
	cout<<"rBn: "<<rBn<<"   I1: "<<I1;
	misc1<<rBn<<"  "<<real(I1);
	I2=0.;
	for (n2=0;n2<dim2->num_puntos;n2++) {
		theta=(dim2->a)+((dim2->b)-(dim2->a))*((dim2->puntos[n2])+1.)/2.;
		seno=sin(theta);
		coseno=cos(theta);
		rdx=0.5*(rAp*seno);
		rdy=0.;
		rdz=0.5*(rAp*coseno+km*rBn);
		rd=sqrt(rdx*rdx+rdy*rdy+rdz*rdz);
		rpnx=-rAp*seno;
		rpny=0.;
		rpnz=rBn*km-rAp*coseno;
		rpn=sqrt(rpnx*rpnx+rpny*rpny+rpnz*rpnz);
		coseno_d=rdz/rd;
		phi_d=atan2(rdy,rdx);
		gl=interpola_cmpx(g->wf,g->r,rd,g->puntos);
		ud=interpola_cmpx(u->wf,u->r,rpn,u->puntos);
		coupling=FuncionAngular3(lp,ld,l,0,coseno,coseno_d,phi,0.);
//		cout<<FuncionAngular3(lp,ld,l,m,coseno,coseno_d,phi,0.)<<"  "<<SphericalHarmonic(l,-m,1.,0.)<<endl;
		//					I1+=seno*senoBn*fl*gl*ud*vpn*coupling*
		//							dim2->pesos[n2]*dim2->pesos[n3]*dim2->pesos[n4]*dim2->pesos[n5]/(rd);
		I2+=seno*fl*gl*ud*vpn*coupling*dim2->pesos[n2]/(rd);
	}
	I2=I2*2.*pow(PI,2.5)/(sqrt(2*l+1.));
	cout<<"    I2: "<<I2<<endl;
	misc1<<"     "<<-real(I2)<<endl;
}
void SourceZR(complejo* rho,distorted_wave* f,distorted_wave* g,int l,double rBn,parametros* parm)
{
	int m,lp,ld;
	double alpha;
	complejo fl,gl,ud;
	lp=g->l;
	ld=f->l;
	alpha=(parm->m_A+1.)/parm->m_A;
	fl=interpola_cmpx(f->wf,f->r,rBn,f->puntos);
	gl=interpola_cmpx(g->wf,g->r,alpha*rBn,g->puntos);
//	if(lp==10) misc1<<rBn<<"  "<<real(fl)/(parm->k_Aa*rBn)<<"  "<<real(gl)/(parm->k_Bb*rBn)<<endl;
	for(m=0;m<=l;m++){
		rho[m]=ClebsGordan(lp,-m,ld,0,l,-m)*ClebsGordan(lp,0,ld,0,l,0)*fl*gl/(alpha*rBn*alpha*rBn);
//		misc1<<rBn<<"  "<<abs(rho[m])<<"  "<<ClebsGordan(lp,-m,ld,0,l,-m)<<"  "<<
//				ClebsGordan(lp,0,ld,0,l,0)<<"  "<<abs(fl)<<"  "<<abs(gl)<<"  "<<abs(alpha*rBn*alpha*rBn)<<endl;
	}
}
void SourceIntegrand(distorted_wave* f,distorted_wave* g,estado* u,potencial_optico* v,potencial_optico* optico,
		 potencial_optico* core,int l,double rBn,parametros* parm)
{
	int n1,n2,m,lp,ld;
	double rAp,theta,rd,rdx,rdz,rpnx,rpnz,rpn,seno,coseno,coseno_d,
	          km,rAn,step;
	complejo fl,gl,ud,coupling,corepot,inpot,remnant,vpn;
	complejo* suma=new complejo[l+1];
	complejo* sumanon=new complejo[l+1];
	ld=f->l;
	lp=g->l;
	km=(parm->m_A+1.)/parm->m_A;
	for(m=0;m<=l;m++){
		suma[m]=0.;
		sumanon[m]=0.;
	}
	rAn=km*rBn;
	vpn=interpola_cmpx(v->pot,v->r,rAn,v->puntos);
	remnant=0.;
	theta=0.;
	step=parm->radio/double(parm->puntos);
	for (n1 =0;n1<parm->puntos; n1++) {
		rAp=(n1+1)*step;
		fl=interpola_cmpx(f->wf,f->r,rAp,f->puntos);
		if(parm->remnant==1) corepot=interpola_cmpx(core->pot,core->r,rAp,core->puntos);
		seno=sin(theta);
		coseno=cos(theta);
		rdx=0.5*(rAp*seno);
		rdz=0.5*(rAp*coseno+km*rBn);
		rd=sqrt(rdx*rdx+rdz*rdz);
		if(parm->remnant==1) inpot=interpola_cmpx(optico->pot,optico->r,rd,optico->puntos);
		rpnx=-rAp*seno;
		rpnz=rBn*km-rAp*coseno;
		rpn=sqrt(rpnx*rpnx+rpnz*rpnz);
		coseno_d=rdz/rd;
		gl=interpola_cmpx(g->wf,g->r,rd,g->puntos);
		ud=interpola_cmpx(u->wf,u->r,rpn,u->puntos);
		coupling=FuncionAngular2(lp,ld,l,coseno,coseno_d);
		if (parm->remnant==1) remnant=inpot-corepot;
		suma[0]=rAp*fl*gl*ud*(vpn-remnant)*coupling/(rd);
		sumanon[0]=rAp*fl*gl*ud*coupling/(rd);
		misc3<<rAp<<"  "<<real(suma[0])<<"  "<<imag(suma[0])<<"  "<<abs(suma[0])<<endl;
	}
	delete[] suma;
	delete[] sumanon;
}
complejo GreenIntegrando(int pts,complejo**** rho,distorted_wave* fl,distorted_wave* Pl,
		parametros* parm,double rBn,int l,int lp,int ld)
{
	int n,m;
	double rBnp,step;
	complejo fl_rBn,Pl_rBn,fl_rBnp,Pl_rBnp;
	fl_rBn=interpola_cmpx(fl->wf,fl->r,rBn,fl->puntos);
	Pl_rBn=interpola_cmpx(Pl->wf,Pl->r,rBn,Pl->puntos);
	complejo suma;
	step=double(parm->radio/pts);
	suma=0.;
	for (n =0;n<pts; n++) {
		rBnp=step*(n+1);
		if(rBn>rBnp) {
			fl_rBnp=interpola_cmpx(fl->wf,fl->r,rBnp,fl->puntos);
			for(m=0;m<=l;m++)
			{
				suma+=real(Pl_rBn*rho[n][l][0][lp]*fl_rBnp*rBnp)*step;
				misc3<<rBnp<<"  "<<real(suma)<<"   "<<endl;
			}
		}
		if(rBn<=rBnp) {
			Pl_rBnp=interpola_cmpx(Pl->wf,Pl->r,rBnp,Pl->puntos);
			for(m=0;m<=l;m++)
			{
				suma+=real(fl_rBn*rho[n][l][0][lp]*Pl_rBnp*rBnp)*step;
				misc3<<rBnp<<"  "<<real(suma)<<"   "<<endl;
			}
		}
	}
	return suma;
}
complejo NeutronWave(complejo* phi,complejo**** rho,distorted_wave* fl,distorted_wave* Pl,
		parametros_integral* dim,parametros* parm,double rBn,int l,int lp,int ld,complejo wronskiano)
{
	int n,m;
	double rBnp;
	complejo* suma=new complejo[l+1];
	complejo* suma2=new complejo[l+1];
	complejo fl_rBn,Pl_rBn,fl_rBnp,Pl_rBnp,rho_int;
	fl_rBn=interpola_cmpx(fl->wf,fl->r,rBn,fl->puntos);
	Pl_rBn=interpola_cmpx(Pl->wf,Pl->r,rBn,Pl->puntos);
	for(m=0;m<=l;m++)
	{
		suma[m]=0.;
		suma2[m]=0.;
	}
	for (n =0;n<dim->num_puntos; n++) {
		rBnp=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[n])+1.)/2.;
		if(rBn>rBnp) {
			fl_rBnp=interpola_cmpx(fl->wf,fl->r,rBnp,fl->puntos);
			for(m=0;m<=l;m++)
			{
				suma[m]+=Pl_rBn*rho[n][l][m][lp]*fl_rBnp*rBnp*dim->pesos[n];
			}
			//			if(l==0 && lp==10) misc3<<rBnp<<"  "<<abs(rho[n][l][0][10])<<endl;
		}
		if(rBn<=rBnp) {
			Pl_rBnp=interpola_cmpx(Pl->wf,Pl->r,rBnp,Pl->puntos);
			for(m=0;m<=l;m++)
			{
				suma2[m]+=fl_rBn*rho[n][l][m][lp]*Pl_rBnp*rBnp*dim->pesos[n];
			}
		}
	}
	for(m=0;m<=l;m++)
	{
		phi[m]=(suma[m]+suma2[m])*((dim->b)-(dim->a))*(1./(wronskiano))/2.;
	}
	delete[] suma;
	delete[] suma2;
}
void NeutronWaveGF(complejo* phi,complejo**** rho,complejo** green,double* r,int puntos_r,
		parametros_integral* dim,parametros* parm,double rBn,int l,int lp,int ld,complejo wronskiano)
{
  int n,m;
  double rBnp;
  complejo* suma=new complejo[l+1];
  complejo greenint,rho_int;
  for(m=0;m<=l;m++)
    {
      suma[m]=0.;
    }
  for (n =0;n<dim->num_puntos; n++) {
    rBnp=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[n])+1.)/2.;
    greenint=interpola2D_cmpx(green,r,r,rBn,rBnp,puntos_r,puntos_r);
    for(m=0;m<=l;m++)
      {
	suma[m]+=greenint*rho[n][l][m][lp]*rBnp*dim->pesos[n];
      }
  }
  for(m=0;m<=l;m++)
    {
      phi[m]=(suma[m])*((dim->b)-(dim->a))/2.;
    }
  delete[] suma;
}
/*****************************************************************************
Interpolacion lineal de una funcion real de 2 variables
 *****************************************************************************/
double interpola2D_dbl(double** funcion,double* r1,double* r2,
		double posicion1,double posicion2,int puntos1,int puntos2)
{
	int indice1,indice2;
	double f11, f12, f21,f22,delta_r1,delta_r2;
	if ((puntos1<3)||(puntos2<3)) Error("El numero de puntos debe ser mayor que 3!");
	delta_r1=r1[puntos1-1]-r1[puntos1-2];
	indice1 = int(ceil(posicion1/delta_r1)) - 1;
	delta_r2=r2[puntos2-1]-r2[puntos2-2];
	indice2 = int(ceil(posicion2/delta_r2)) - 1;
	if (indice1 > puntos1 - 2)
		indice1=puntos1-2;
	if (indice1 <= 0)
		indice1=1;
	if (indice2 > puntos2 - 2)
		indice2=puntos2-2;
	if (indice2 <= 0)
		indice2=1;
	f11=funcion[indice1][indice2];
	f12=funcion[indice1][indice2+1];
	f21=funcion[indice1+1][indice2];
	f22=funcion[indice1+1][indice2+1];
	return (f11*(r1[indice1+1]-posicion1)*(r2[indice2+1]-posicion2)+f21*(posicion1-r1[indice1])*(r2[indice2+1]-posicion2)
			+f12*(r1[indice1+1]-posicion1)*(posicion2-r2[indice2])+f22*(posicion1-r1[indice1])*(posicion2-r2[indice2]))/(delta_r1*delta_r2);
}
/*****************************************************************************
Interpolacion lineal de una funcion compleja de 2 variables
 *****************************************************************************/
complejo interpola2D_cmpx(complejo** funcion,double* r1,double* r2,
		double posicion1,double posicion2,int puntos1,int puntos2)
{
	int indice1,indice2;
	complejo f11, f12, f21,f22;
	double delta_r1,delta_r2;
	if ((puntos1<3)||(puntos2<3)) Error("El numero de puntos debe ser mayor que 3!");
	delta_r1=r1[puntos1-1]-r1[puntos1-2];
	indice1 = int(ceil(posicion1/delta_r1)) - 1;
	delta_r2=r2[puntos2-1]-r2[puntos2-2];
	indice2 = int(ceil(posicion2/delta_r2)) - 1;
	if (indice1 > puntos1 - 2)
		indice1=puntos1-2;
	if (indice1 <= 0)
		indice1=1;
	if (indice2 > puntos2 - 2)
		indice2=puntos2-2;
	if (indice2 <= 0)
		indice2=1;
	f11=funcion[indice1][indice2];
	f12=funcion[indice1][indice2+1];
	f21=funcion[indice1+1][indice2];
	f22=funcion[indice1+1][indice2+1];
//	misc2<<posicion1<<"  "<<posicion2<<"  "<<delta_r1<<"  "<<delta_r2<<"  "<<indice1<<"  "<<indice2<<"  "<<f11<<"  "<<f21<<"  "<<f12<<"  "<<f22<<"  "<<endl;
//	exit(0);
	return (f11*(r1[indice1+1]-posicion1)*(r2[indice2+1]-posicion2)+f21*(posicion1-r1[indice1])*(r2[indice2+1]-posicion2)
			+f12*(r1[indice1+1]-posicion1)*(posicion2-r2[indice2])+f22*(posicion1-r1[indice1])*(posicion2-r2[indice2]))/(delta_r1*delta_r2);
}
complejo NeutronWaveTest(complejo* phi,complejo* rho,distorted_wave* fl,distorted_wave* Pl,
		parametros_integral* dim,parametros* parm,double rBn,int l,int lp,int ld,complejo wronskiano)
{
	int n,m;
	double rBnp;
	complejo suma;
	complejo suma2;
	complejo fl_rBn,Pl_rBn,fl_rBnp,Pl_rBnp,rho_int;
	fl_rBn=interpola_cmpx(fl->wf,fl->r,rBn,fl->puntos);
	Pl_rBn=interpola_cmpx(Pl->wf,Pl->r,rBn,Pl->puntos);
	suma=0.;
	suma2=0.;
//	cout<<"rBn: "<<rBn<<endl;
	for (n =0;n<dim->num_puntos; n++) {
		rBnp=(dim->a)+(rBn-(dim->a))*((dim->puntos[n])+1.)/2.;
		rho_int=interpola_cmpx(rho,fl->r,rBnp,fl->puntos);
		fl_rBnp=interpola_cmpx(fl->wf,fl->r,rBnp,fl->puntos);
		suma+=Pl_rBn*rho_int*fl_rBnp*rBnp*dim->pesos[n];
//		misc3<<rBnp<<"  "<<real(Pl_rBn*fl_rBnp/(wronskiano))<<"  "<<imag(Pl_rBn*fl_rBnp/(wronskiano))<<endl;
	}
	for (n =0;n<dim->num_puntos; n++) {
		rBnp=(rBn)+((dim->b)-rBn)*((dim->puntos[n])+1.)/2.;
		rho_int=interpola_cmpx(rho,fl->r,rBnp,fl->puntos);
		Pl_rBnp=interpola_cmpx(Pl->wf,Pl->r,rBnp,Pl->puntos);
		suma2+=fl_rBn*rho_int*Pl_rBnp*rBnp*dim->pesos[n];
//		misc3<<rBnp<<"  "<<real(Pl_rBnp*fl_rBn/(wronskiano))<<"  "<<imag(Pl_rBnp*fl_rBn/(wronskiano))<<endl;
	}
	phi[0]=(suma*(rBn-(dim->a))*0.5+suma2*((dim->b)-rBn)*0.5)*(1./(wronskiano));
}
complejo NeutronWaveResonant(complejo* phi,complejo**** rho, estado* st,double En, double E, double absorcion,
		parametros_integral* dim,parametros* parm,double rBn,int l,int lp)
{
	int n,m;
	double rBnp;
	complejo lorentz;
	complejo* suma=new complejo[l+1];
	complejo stp,stint;
	stp=interpola_cmpx(st->wf,st->r,rBn,st->puntos);
	lorentz=(HC*HC/(2.*AMU))*(1./((En-E)+I*absorcion));
	for(m=0;m<=l;m++)
	{
		suma[m]=0.;
	}
	for (n =0;n<dim->num_puntos; n++) {
		rBnp=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[n])+1.)/2.;
		stint=interpola_cmpx(st->wf,st->r,rBnp,st->puntos);
		for(m=0;m<=l;m++)
		{
			suma[m]+=rho[n][l][m][lp]*stint*rBnp*rBnp*dim->pesos[n];
		}
	}
	for(m=0;m<=l;m++)
	{
		phi[m]=(suma[m]*((dim->b)-(dim->a))*0.5)*stp*rBn*lorentz;
	}
	delete[] suma;
}
void ElasticBreakup(complejo*** T,complejo**** rho,double En,potencial_optico* v_up,potencial_optico* v_down,
		parametros_integral* dim,parametros* parm,int l,int lp,double kn)
{
	int n,m;
		double rBnp,redfac,masaT,masa_res,carga_trans,masa_trans,up,down,spin;
		redfac=2.*AMU/(HC*HC);
		complejo* suma=new complejo[l+1];
		complejo stintup,stintdown;
		distorted_wave* chi_lup=new distorted_wave;
		distorted_wave* chi_ldown=new distorted_wave;
		spin=0.;
		ofstream fp("neutron_distorted_wave.txt");
		chi_lup->energia=En;
		chi_lup->l=l;
		chi_lup->spin=spin;
		chi_lup->j=l+spin;
		chi_ldown->energia=En;
		chi_ldown->l=l;
		chi_ldown->spin=spin;
		chi_ldown->j=l-spin;
		masaT=parm->T_carga+parm->T_N;
		masa_res=parm->res_carga+parm->res_N;
		carga_trans=parm->res_carga-parm->T_carga;
		masa_trans=masa_res-masaT;
		up=((l+1.)/sqrt((l+1.)*(l+1.)+l*l));
		down=(l/sqrt((l+1.)*(l+1.)+l*l));
		//fp<<endl<<"& spin up"<<endl;
		GeneraDWspin(chi_lup,v_up,0.,masa_trans*masaT/(masa_trans+masaT),
				parm->radio,parm->puntos,parm->matching_radio,&fp);
		//fp<<"& spin down"<<endl;
		GeneraDWspin(chi_ldown,v_down,0.,masa_trans*masaT/(masa_trans+masaT),
				parm->radio,parm->puntos,parm->matching_radio,&fp);
		for(m=0;m<=l;m++)
		{
			suma[m]=0.;
		}
		for (n =0;n<dim->num_puntos; n++) {
			rBnp=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[n])+1.)/2.;
			stintup=up*interpola_cmpx(chi_lup->wf,chi_lup->r,rBnp,chi_lup->puntos);
			stintdown=down*interpola_cmpx(chi_ldown->wf,chi_ldown->r,rBnp,chi_ldown->puntos);
			for(m=0;m<=l;m++)
			{
				suma[m]+=rho[n][l][m][lp]*(stintup+stintdown)*rBnp*dim->pesos[n];
			}
//			if (lp==0) misc2<<rBnp<<"  "<<abs(stintup)<<"  "<<abs(stintdown)<<"  "<<abs(suma[0])<<endl;
		}
		for(m=0;m<=l;m++)
		{
			T[l][m][lp]=4.*PI*(suma[m]*((dim->b)-(dim->a))*0.5)/(kn*redfac);
		}
		delete[] suma;
		delete[] chi_lup;
		delete[] chi_ldown;
}
complejo NeutronWaveResonantTest(complejo* phi,complejo* rho, estado* st,double En, double E, double absorcion,
		parametros_integral* dim,parametros* parm,double rBn)
{
	int n,m;
	double rBnp;
	complejo lorentz;
	complejo suma;
	complejo stp,stint,fl_rBnp,Pl_rBnp,rho_int;
	stp=interpola_cmpx(st->wf,st->r,rBn,st->puntos);
	lorentz=(HC*HC/(2.*AMU))*(1./((En-E)+I*absorcion));
//	cout<<(HC*HC/(2.*AMU))<<"  absorcion: "<<absorcion<<"  cociente: "<<(1./((En-E)+I*absorcion))
//		<<"  diferencia: "<<(En-E)<<"  En: "<<En<<"  E: "<<E<<"  lorentz: "<<lorentz<<endl;
	suma=0.;
//	cout<<"rBn: "<<rBn<<endl;
	for (n =0;n<dim->num_puntos; n++) {
		rBnp=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[n])+1.)/2.;
		rho_int=interpola_cmpx(rho,st->r,rBnp,st->puntos);
		stint=interpola_cmpx(st->wf,st->r,rBnp,st->puntos);
		suma+=rho_int*stint*rBnp*dim->pesos[n];
//		misc3<<rBnp<<"  "<<real(Pl_rBn*fl_rBnp/(wronskiano))<<"  "<<imag(Pl_rBn*fl_rBnp/(wronskiano))<<endl;
	}
	phi[0]=(suma*((dim->b)-(dim->a))*0.5)*stp*lorentz;
}
double Absorcion(potencial_optico* pot,complejo**** wf,parametros_integral* dim,int l,int lmax)
{
	int n,m,lp;
	double R;
	double sumam=0.;
	complejo pot_int;
	for(lp=0;lp<lmax;lp++)
	{
		for(n=0;n<dim->num_puntos;n++)
		{
			R=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[n])+1.)/2.;
			pot_int=interpola_cmpx(pot->pot,pot->r,R,pot->puntos);
			for(m=0;m<=l;m++)
			{
				sumam+=-R*R*imag(pot_int)*abs(wf[n][l][m][lp])*abs(wf[n][l][m][lp])*dim->pesos[n]*((dim->b)-(dim->a))/2.;
			}
//			if(n==10) misc3<<l<<"  "<<lp<<"  "<<abs(suma[0])<<"  "<<abs(abs(wf[n][l][0][lp]))<<endl;
		}
//		misc1<<lp<<"  "<<abs(sumam)<<endl;
	}
	return sumam;
}
double Absorcion2(potencial_optico* pot,estado* wf)
{
	int n,m,lp;
	double R;
	double suma=0.;
	parametros_integral *dim=new parametros_integral;
	complejo pot_int,st;
	dim->a=0.;
	dim->b=30.;
	dim->num_puntos=50;
	GaussLegendre(dim->puntos,dim->pesos,dim->num_puntos);
	for(n=0;n<dim->num_puntos;n++)
	{
		R=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[n])+1.)/2.;
		pot_int=interpola_cmpx(pot->pot,pot->r,R,pot->puntos);
		st=interpola_cmpx(wf->wf,wf->r,R,wf->puntos);
		suma+=-R*R*imag(pot_int)*abs(st)*abs(st)*dim->pesos[n]*((dim->b)-(dim->a))/2.;
//		cout<<R<<"  "<<suma<<"  "<<pot_int<<"  "<<st<<endl;
	}
	delete[] dim;
	return suma;
}
double AbsorcionPrior(double* direct,double* non_orth,double* cross,
		potencial_optico* pot,complejo**** wf,complejo**** non,parametros_integral* dim,int l,int lmax,double r_F)
{
	int n,m,lp;
	double R,suma,sumaUT,sumaHM;
	complejo pot_int,UT,HM;
	suma=0.;
    direct[0]=0.;
    non_orth[0]=0.;
    cross[0]=0.;
    suma=0.;
    sumaUT=0.;
    sumaHM=0.;
//    misc2<<"Quillo!!+++++++++++++++++++++++++++++++++++++++++++"<<endl;
    for(lp=0;lp<lmax;lp++)
    {
    	for(m=0;m<=lp;m++)
    	{
    		for(n=0;n<dim->num_puntos;n++)
    		{
    			R=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[n])+1.)/2.;
    			pot_int=interpola_cmpx(pot->pot,pot->r,R,pot->puntos);
    			if(m==0){
    				UT=wf[n][l][0][lp];
    				HM=non[n][l][0][lp];
    			}
    			if(m>0){
    				UT=sqrt(2.)*wf[n][l][m][lp];
    				HM=sqrt(2.)*non[n][l][m][lp];
    			}
    			if (R<r_F)
    			{
    				suma+=-imag(pot_int)*abs(UT-HM)*abs(UT-HM)*dim->pesos[n]*((dim->b)-(dim->a))/2.;
    				sumaUT+=-imag(pot_int)*abs(UT)*abs(UT)*dim->pesos[n]*((dim->b)-(dim->a))/2.;
    				sumaHM+=-imag(pot_int)*abs(HM)*abs(HM)*dim->pesos[n]*((dim->b)-(dim->a))/2.;
    			}
    		}

    	}
    }
    return suma;
}
double AbsorcionNL(complejo** pot,complejo**** wf,complejo**** non,parametros_integral* dim,int l,int lmax
		   ,double* r,int puntos_r)
{
  int n,nn,m,lp;
  double R,RR;
  complejo pot_int,UT,HM,UTT,HMM,suma,sumaUT,sumaHM;
  suma=0.;
  sumaUT=0.;
  sumaHM=0.;
  for(lp=0;lp<lmax;lp++)
    {
      for(m=0;m<=lp;m++)
    	{
	  for(n=0;n<dim->num_puntos;n++)
	    {
	      R=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[n])+1.)/2.;
	      if(m==0){
		UT=wf[n][l][0][lp];
		HM=non[n][l][0][lp];
	      }
	      if(m>0){
		UT=sqrt(2.)*wf[n][l][m][lp];
		HM=sqrt(2.)*non[n][l][m][lp];
	      }
	      for(nn=0;nn<dim->num_puntos;nn++)
		{
		  RR=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[nn])+1.)/2.;
		  pot_int=interpola2D_cmpx(pot,r,r,R,RR,puntos_r,puntos_r);
		  if(m==0){
		    UTT=wf[nn][l][0][lp];
		    HMM=non[nn][l][0][lp];
		  }
		  if(m>0){
		    UTT=sqrt(2.)*wf[nn][l][m][lp];
		    HMM=sqrt(2.)*non[nn][l][m][lp];
		  }
		  suma+=-imag(pot_int)*(UT-HM)*conj(UTT-HMM)*dim->pesos[n]*((dim->b)-(dim->a))*dim->pesos[nn]*((dim->b)-(dim->a))/4.;
		  sumaUT+=-imag(pot_int)*(UT)*conj(UT)*dim->pesos[n]*((dim->b)-(dim->a))*dim->pesos[nn]*((dim->b)-(dim->a))/4.;
		  sumaHM+=-imag(pot_int)*(HM)*conj(HM)*dim->pesos[n]*((dim->b)-(dim->a))*dim->pesos[nn]*((dim->b)-(dim->a))/4.;
		}
	    }
    	}
    }
  return abs(suma);
}

double ElasticBreakupCross(complejo*** Teb,int l,int lmax)
{
	int m,lp;
	double suma;
	suma=0.;
    for(lp=0;lp<lmax;lp++)
    {
    	suma+=abs(Teb[l][0][lp])*abs(Teb[l][0][lp]);
    	for(m=1;m<=lp;m++)
    	{
    		suma+=2.*abs(Teb[l][m][lp])*abs(Teb[l][m][lp]);
    	}
    }
	return suma;
}
double AbsorcionPriorTest(double* direct,double* non_orth,double* cross,double* r,int puntos,
		potencial_optico* pot,complejo* wf,complejo* non,parametros_integral* dim)
{
	int n,m,lp;
	double R,suma;
	complejo pot_int,wf_int,non_int;
	suma=0.;
	direct[0]=0.;
	non_orth[0]=0.;
	cross[0]=0.;
	for(n=0;n<dim->num_puntos;n++)
	{
		R=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[n])+1.)/2.;
		pot_int=interpola_cmpx(pot->pot,pot->r,R,pot->puntos);
		wf_int=interpola_cmpx(wf,r,R,puntos);
		non_int=interpola_cmpx(non,r,R,puntos);
		direct[0]+=-imag(pot_int)*abs(wf_int)*abs(wf_int)*dim->pesos[n]*((dim->b)-(dim->a))/2.;
		non_orth[0]+=-imag(pot_int)*abs(non_int)*abs(non_int)*dim->pesos[n]*((dim->b)-(dim->a))/2.;
		cross[0]+=-2.*imag(pot_int)*real(wf_int*non_int)*dim->pesos[n]*((dim->b)-(dim->a))/2.;
		suma+=direct[0]+non_orth[0]-cross[0];
	}
	return direct[0]+non_orth[0]-cross[0];
}
double AbsorcionAngular(potencial_optico* pot,complejo**** wf,complejo**** non,parametros_integral* dim,parametros* parm,
		double theta, double* direct, double* non_orth, double* cross, double* cross_j)
{
	int n,m,lp,l;
	double R,suma,armonico,costheta;
	complejo pot_int,UT,HM,phase2;
	costheta=cos(theta);
	suma=0.;
	for(n=0;n<dim->num_puntos;n++)
	{
		R=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[n])+1.)/2.;
		pot_int=interpola_cmpx(pot->pot,pot->r,R,pot->puntos);
		for(l=0;l<parm->ltransfer;l++)
		{
			UT=0.;
			HM=0.;
			for(lp=0;lp<parm->lmax;lp++)
			{
				armonico=gsl_sf_legendre_sphPlm(lp,0,costheta);
				UT+=wf[n][l][0][lp]*armonico;
				HM+=non[n][l][0][lp]*armonico;
			}
			suma+=-imag(pot_int)*abs(UT-HM)*abs(UT-HM)*dim->pesos[n]*((dim->b)-(dim->a))/2.;
//			suma+=-imag(pot_int)*abs(UT)*abs(UT)*dim->pesos[n]*((dim->b)-(dim->a))/2.;
			for(m=1;m<parm->lmax;m++)
			{
				UT=0.;
				HM=0.;
				for(lp=m;lp<parm->lmax;lp++)
				{
					armonico=gsl_sf_legendre_sphPlm(lp,m,costheta);
					UT+=wf[n][l][m][lp]*armonico;
					HM+=non[n][l][m][lp]*armonico;
				}
				suma+=-2.*imag(pot_int)*abs(UT-HM)*abs(UT-HM)*dim->pesos[n]*((dim->b)-(dim->a))/2.;
//				suma+=-2.*imag(pot_int)*abs(UT)*abs(UT)*dim->pesos[n]*((dim->b)-(dim->a))/2.;
			}
		}
	}
	return suma;
}

double ElasticBreakupAngular(complejo*** Teb,int lmax,double theta)
{
	int m,lp,l;
	double armonico,costheta,B;
	complejo A;
	costheta=cos(theta);
	B=0.;
	for(l=0;l<lmax;l++)
	{
		A=0.;
		for(lp=0;lp<lmax;lp++)
		{
			armonico=gsl_sf_legendre_sphPlm(lp,0,costheta);
			A+=Teb[l][0][lp]*armonico;
		}
		B+=abs(A)*abs(A);
		for(m=1;m<=lmax;m++)
		{
			A=0.;
			for(lp=m;lp<lmax;lp++)
			{
				armonico=gsl_sf_legendre_sphPlm(lp,m,costheta);
				A+=Teb[l][m][lp]*armonico;
			}
			B+=2.*abs(A)*abs(A);
		}
	}
	return B;
}
double TotalBreakup(complejo**** wf,complejo**** rho,parametros* parm, parametros_integral* dim,int l)
{
	int m,lp,n;
	double R;
	double suma;
	suma=0.;
    for(lp=0;lp<parm->lmax;lp++)
    {
    	for(m=0;m<=lp;m++)
    	{
    		for(n=0;n<dim->num_puntos;n++)
    		{
    			R=(dim->a)+((dim->b)-(dim->a))*((dim->puntos[n])+1.)/2.;
    				suma+=imag(wf[n][l][m][lp]*conj(rho[n][l][m][lp]))*R*R*dim->pesos[n]*((dim->b)-(dim->a))/2.;
    				if(m>0) suma+=imag(wf[n][l][m][lp]*conj(rho[n][l][m][lp]))*R*R*dim->pesos[n]*((dim->b)-(dim->a))/2.;
    		}
    	}
    }
	return suma;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                       //
//     Funcion <lp -m ld 0 | l -m> [Y^lp(costheta)Y^ld(costheta_d)]^K_M Y^l_-m(costheta_Bn)              //
//                                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////////////////////
complejo FuncionAngular(int lp,int ld,int l,int m,int K,int M,double costheta, double costheta_d, double costheta_Bn)
{
	if (l>lp+ld) Error("l>lp+ld en FuncionAngular");
	if (l<abs(lp-ld)) Error("l<lp-ld en FuncionAngular");
	if (m>l) Error("m>l en FuncionAngular");
	if (m>lp) Error("m>lp en FuncionAngular");
	int mlp;
	complejo suma=0.;
	if ((M>=0)&&(abs(M<=ld))) suma+=ClebsGordan(lp,0,ld,M,K,M)*gsl_sf_legendre_sphPlm(lp,0,costheta)*
			gsl_sf_legendre_sphPlm(ld,M,costheta_d);
	if ((M<0)&&(abs(M<=ld))) suma+=pow(-1.,M)*ClebsGordan(lp,0,ld,M,K,M)*gsl_sf_legendre_sphPlm(lp,0,costheta)*
			gsl_sf_legendre_sphPlm(ld,M,costheta_d);
	for (mlp=1;mlp<=lp;mlp++)
	{
		if ((M-mlp>=0)&&(abs(M-mlp)<=ld)) suma+=ClebsGordan(lp,mlp,ld,M-mlp,K,M)*gsl_sf_legendre_sphPlm(lp,mlp,costheta)*
				gsl_sf_legendre_sphPlm(ld,M-mlp,costheta_d);
		if ((M-mlp<0)&&(abs(M-mlp)<=ld)) suma+=pow(-1.,M-mlp)*ClebsGordan(lp,mlp,ld,M-mlp,K,m)*gsl_sf_legendre_sphPlm(lp,mlp,costheta)*
				gsl_sf_legendre_sphPlm(ld,mlp-M,costheta_d);
		if ((M+mlp)>=0&&(abs(M+mlp)<=ld)) suma+=pow(-1.,mlp)*ClebsGordan(lp,-mlp,ld,m+mlp,K,M)*gsl_sf_legendre_sphPlm(lp,mlp,costheta)*
				gsl_sf_legendre_sphPlm(ld,M+mlp,costheta_d);
		if ((M+mlp<0)&&(abs(M+mlp)<=ld)) suma+=pow(-1.,M)*ClebsGordan(lp,-mlp,ld,m+mlp,K,M)*gsl_sf_legendre_sphPlm(lp,mlp,costheta)*
				gsl_sf_legendre_sphPlm(ld,-mlp-M,costheta_d);
	}
//suma=1.;
	if (-m>=0) suma*=gsl_sf_legendre_sphPlm(l,-m,costheta_Bn);
	if (-m<0) suma*=pow(-1.,m)*gsl_sf_legendre_sphPlm(l,m,costheta_Bn);
	return ClebsGordan(lp,-m,ld,0,l,-m)*suma;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                       //
//     Funcion [Y^lp(costheta)Y^ld(costheta_d)]^l_0                                                      //
//                                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////////////////////
complejo FuncionAngular2(int lp,int ld,int l,double costheta, double costheta_d)
{
	if (l>lp+ld) Error("l>lp+ld en FuncionAngular2");
	if (l<abs(lp-ld)) Error("l<lp-ld en FuncionAngular2");
	int mlp;
	complejo suma=0.;
	suma+=ClebsGordan(lp,0,ld,0,l,0)*gsl_sf_legendre_sphPlm(lp,0,costheta)*
			gsl_sf_legendre_sphPlm(ld,0,costheta_d);
	for (mlp=1;mlp<=lp;mlp++)
	{
		if(mlp<=ld) suma+=pow(-1.,mlp)*gsl_sf_legendre_sphPlm(lp,mlp,costheta)*
				gsl_sf_legendre_sphPlm(ld,mlp,costheta_d)*(ClebsGordan(lp,mlp,ld,-mlp,l,0)+ClebsGordan(lp,-mlp,ld,mlp,l,0));
	}

	return suma;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                       //
//     Funcion [Y^lp(costheta,phi)Y^ld(costheta_d,phi_d)]^l_m                                            //
//                                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////////////////////
complejo FuncionAngular3(int lp,int ld,int l,int m,double costheta, double costheta_d,double phi, double phi_d)
{
	if (l>lp+ld) Error("l>lp+ld en FuncionAngular3");
	if (l<abs(lp-ld)) Error("l<lp-ld en FuncionAngular3");
	int mlp;
	complejo suma=0.;
	if(abs(m)<=ld) suma+=ClebsGordan(lp,0,ld,m,l,m)*SphericalHarmonic(lp,0,costheta,phi)*SphericalHarmonic(ld,m,costheta_d,phi_d);
	for (mlp=1;mlp<=lp;mlp++)
	{
		if(abs(m-mlp)<=ld) suma+=ClebsGordan(lp,mlp,ld,m-mlp,l,m)*SphericalHarmonic(lp,mlp,costheta,phi)*
				SphericalHarmonic(ld,m-mlp,costheta_d,phi_d);
		if(abs(m+mlp)<=ld) suma+=ClebsGordan(lp,-mlp,ld,m+mlp,l,m)*SphericalHarmonic(lp,-mlp,costheta,phi)*
								SphericalHarmonic(ld,m+mlp,costheta_d,phi_d);
	}
	return suma;
}
complejo SphericalHarmonic(int l, int m, double costheta, double phi)
{
	complejo valor;
	if (m>l) {cout<<"Cuidado, m>l en SphericalHarmonic. Se devuelve 0."<<endl; return 0.;}
	if (m>=0) valor=gsl_sf_legendre_sphPlm(l,m,costheta)*exp(I*double(m)*phi);
	if (m<0) valor=pow(-1.,m)*gsl_sf_legendre_sphPlm(l,-m,costheta)*exp(I*double(m)*phi);
	return valor;
}
complejo GeneraGreenFunctionLigada(distorted_wave *regular,distorted_wave *irregular,potencial_optico *potencial,
		double radio_max,int puntos,double q1q2,double masa,double spin) {
	int ND,i,wronsk_point;
	complejo *sx=new complejo[puntos],*vs=new complejo[puntos];
	complejo *v=new complejo[puntos],*vv=new complejo[puntos];
	estado* autofuncion=new estado;
	double hbarx,dd,Wlim,centr,ls,delta_r,radio,k,Etrial,Emax,Emin,eigen,wronsk_radio;
	complejo derivada_irregular,derivada_regular,wronskiano;
	delta_r=radio_max/double(puntos);
	wronsk_radio=10.;
	wronsk_point=int(ceil(wronsk_radio/delta_r)-1.);
	eigen=-41.0226066052632;
	hbarx=HC*HC/(2.*AMU*masa);
	dd=delta_r*delta_r/hbarx;
	Wlim=1.e-13;
	centr=(regular->l*(regular->l+1.))*hbarx;
	regular->puntos=puntos;
	regular->radio=radio_max;
	irregular->puntos=puntos;
	irregular->radio=radio_max;
	autofuncion->puntos=puntos;
	autofuncion->radio=radio_max;
	ls=(regular->j*(regular->j+1.)-regular->l*(regular->l+1.)-spin*(spin+1.));
	cout<<"ls en GeneraGreenFunctionLigada: "<<ls<<endl;
	// a�ade los t�rminos Coulomb y spin-�rbita *********************************************************************
	for (i=0;i<puntos;i++) {
		regular->r[i]=delta_r*(i+1);
		irregular->r[i]=delta_r*(i+1);
		if(regular->r[i]>potencial->radio_coul) v[i]=potencial->pot[i]+E_CUADRADO*q1q2/regular->r[i];
		if(potencial->r[i]<=potencial->radio_coul) v[i]=potencial->pot[i]
		               +E_CUADRADO*q1q2*(3.-(regular->r[i]/potencial->radio_coul)*
		            		   (regular->r[i]/potencial->radio_coul))/(2.*potencial->radio_coul);
		vs[i]=-2.*(potencial->Vso)*exp((regular->r[i]-potencial->radioso)/potencial->aso)/
				(regular->r[i]*potencial->aso*(1.+exp((regular->r[i]-potencial->radioso)
						/potencial->aso))*(1.+exp((regular->r[i]-potencial->radioso)/potencial->aso)));
		vv[i]=v[i]+ls*vs[i]+(centr)/(regular->r[i]*regular->r[i]);
//		misc1<<regular->r[i]<<"  "<<real(vv[i])<<"  "<<imag(vv[i])<<endl;
	}
	ND=0;
	Etrial=irregular->energia;
	k=sqrt(-2.*AMU*masa*Etrial)/HC;
	irregular->wf[puntos-1]=1.e-10;
	irregular->wf[puntos-2]=1.e-10/(1.-k*delta_r);
	for (i=puntos-2;i>0;i--) {
		irregular->wf[i-1]= (2.*(1.-0.416666667*dd*(-vv[i]+Etrial))
				*irregular->wf[i]-(1.+0.083333333*dd*(-vv[i+1]+Etrial))
				*irregular->wf[i+1])/(1.+0.083333333*dd*(-vv[i-1]+Etrial));
		if (real(irregular->wf[i+1])*real(irregular->wf[i])<0. ) ND=ND+1;
	}
	regular->wf[0]=pow(delta_r,(regular->l+1));
	regular->wf[1]=pow(2.0*delta_r,(regular->l+1));
	autofuncion->wf[0]=pow(delta_r,(regular->l+1));
	autofuncion->wf[1]=pow(2.0*delta_r,(regular->l+1));
	ND=0;
	for (i=1;i<puntos-1;i++) {
		sx[i]=dd*(vv[i]-Etrial);
		regular->wf[i+1]=(2.0+sx[i])*regular->wf[i]-regular->wf[i-1];
		autofuncion->wf[i+1]=(2.0+dd*(vv[i]-eigen))*autofuncion->wf[i]-autofuncion->wf[i-1];
		if (real(regular->wf[i+1])*real(regular->wf[i])<0. ) ND=ND+1;
	}
	NormalizaD(regular,regular,regular->radio,regular->puntos,'s');
	NormalizaD(irregular,irregular,irregular->radio,irregular->puntos,'s');
	Normaliza(autofuncion,autofuncion,irregular->radio,irregular->puntos,'s');
	delta_r=regular->r[wronsk_point]-regular->r[wronsk_point-1];
	for (i=2;i<puntos-3;i++) {
		if(i==wronsk_point)
		{
			derivada_irregular=(-irregular->wf[i+2]+8.*irregular->wf[i+1]-8.*irregular->wf[i-1]+irregular->wf[i-2])/(12.*delta_r);
			derivada_regular=(-regular->wf[i+2]+8.*regular->wf[i+1]-8.*regular->wf[i-1]+regular->wf[i-2])/(12.*delta_r);
			wronskiano=-derivada_irregular*regular->wf[i]+derivada_regular*irregular->wf[i];
		}
	}
	delete[] vv;
	delete[] v;
	delete[] sx;
	delete[] vs;
	return wronskiano;
}
void HulthenWf(estado *st,double radio_max,int puntos)
{
	int i;
	double delta_r,N,gamma,beta;
	delta_r=radio_max/double(puntos);
	st->puntos=puntos;
	N=0.884872872225158;
	gamma=0.2316;
	beta=1.384968;
	for (i=0;i<puntos;i++) {
		st->r[i]=delta_r*(i+1);
		st->wf[i]=N*(exp(-gamma*st->r[i])-exp(-beta*st->r[i]));
		st->wf[i]=st->wf[i]/st->r[i];
	}
	st->radio=radio_max;
	//  Normalizacion
//	Normaliza(st,st,radio_max,puntos,'s');
}
void TalysInput(double* lmenos,double* lmas,double energia_trans,parametros* parm,ofstream* fp,ofstream* fp2,ofstream* fp3,double s)
{
	int l,parity;
	double Ex,factor,J,Jp,j_menos,j_mas,norma;
	double* parity1=new double[parm->ltransfer+int(ceil(s))+2];
	double* parity2=new double[parm->ltransfer+int(ceil(s))+2];
	Ex=energia_trans-parm->B_Sn;
	*fp3<<"NEnergy= "<<Ex<<"   Entries "<<parm->ltransfer+int(ceil(s))+2<<endl;
	parity=1;
	for(l=0;l<parm->ltransfer+int(ceil(s))+2;l++){
		parity1[l]=0.;
		parity2[l]=0.;
	}
	for(l=1;l<parm->ltransfer;l+=2){
		j_menos=l-0.5;
		if(l==0) j_menos=0.5;
		for(Jp=abs(j_menos-s);Jp<=(j_menos+s);Jp++){
			factor=(2.*Jp+1.)/((j_menos+s)*(j_menos+s+2.)-(abs(j_menos-s)-1.)*(abs(j_menos-s)+1.));
			parity1[int(Jp)]+=factor*lmenos[l];
		}
		j_mas=l+0.5;
		for(Jp=abs(j_mas-s);Jp<=(j_mas+s);Jp++){
			factor=(2.*Jp+1.)/((j_mas+s)*(j_mas+s+2.)-(abs(j_mas-s)-1.)*(abs(j_mas-s)+1.));
			parity1[int(Jp)]+=factor*lmas[l];
		}
	}
	for(l=0;l<parm->ltransfer;l+=2){
		j_menos=l-0.5;
		if(l==0) j_menos=0.5;
		for(Jp=abs(j_menos-s);Jp<=(j_menos+s);Jp++){
			factor=(2.*Jp+1.)/((j_menos+s)*(j_menos+s+2.)-(abs(j_menos-s)-1.)*(abs(j_menos-s)+1.));
			parity2[int(Jp)]+=factor*lmenos[l];
		}
		j_mas=l+0.5;
		for(Jp=abs(j_mas-s);Jp<=(j_mas+s);Jp++){
			factor=(2.*Jp+1.)/((j_mas+s)*(j_mas+s+2.)-(abs(j_mas-s)-1.)*(abs(j_mas-s)+1.));
			parity2[int(Jp)]+=factor*lmas[l];
		}
	}
	norma=0.;
	*fp<<Ex<<"  ";
	for(Jp=0;Jp<parm->ltransfer+int(ceil(s))+2;Jp++){
		*fp<<parity1[int(Jp)]<<"    ";
		norma+=parity1[int(Jp)]+parity2[int(Jp)];
	}
	*fp<<endl;
	*fp2<<Ex<<"  ";
	for(Jp=0;Jp<parm->ltransfer+int(ceil(s))+2;Jp++){
		*fp2<<parity2[int(Jp)]<<"    ";
	}
	*fp2<<endl;

	for(Jp=0;Jp<parm->ltransfer+int(ceil(s))+2;Jp++){
		*fp3<<left<<setw(5)<<Jp<<setw(12)<<left<<parity1[int(Jp)]/norma<<left<<setw(12)<<parity2[int(Jp)]/norma<<endl;
	}
}
/*****************************************************************************
Mensaje de error
 *****************************************************************************/
void Error(const char *text)
{
  time_t rawtime;
  time(&rawtime);
  fprintf(stderr,"---------------------------------------------\n");
  fprintf(stderr,"Error:\n");
  fprintf(stderr,"* %s\n",text);
  fprintf(stdout,"tras %3.1f s. , %s\n",(double) ( clock() / CLOCKS_PER_SEC ),
                                    asctime(localtime(&rawtime)));
  fprintf(stderr,"---------------------------------------------\n\n");

  exit(1);
}
/*****************************************************************************
Lee un n�mero decimal
 *****************************************************************************/
void ReadParD(char *s,const char key[20], int *par)
{
	int l,r;
	l = strlen(key);
	if (!strncmp(s,key,l))
	{
		r = sscanf(s,"%*s %d",par);
		if (r<1) { fprintf(stderr,"ERROR: No puede leerse %s en el fichero de parametros (r=%d)\n",key,r); exit(1); }
		informe<<key<<"="<<*par<<endl;
	}
	fflush(stderr);

}
/*****************************************************************************
Lee un float
 *****************************************************************************/
void ReadParF(char *s,const char key[20],double *par)
{
	int l,r;

	l = strlen(key);
	if (!strncmp(s,key,l))
	{
		r = sscanf(s,"%*s %lf",par);
		if (r<1) { fprintf(stderr,"ERROR: No puede leerse %s en el fichero de parametros\n",key); exit(1); }
		informe<<key<<"="<<*par<<endl;
	}
	fflush(stderr);
}
/*****************************************************************************
Lee una cadena de caracteres
 *****************************************************************************/
void ReadParS(char *s,const char key[20], char *par)
{
	int l,r;

	l = strlen(key);
	if (!strncmp(s,key,l))
	{
		r = sscanf(s,"%*s %s",par);
		if (r<1) { fprintf(stderr,"%s",key); par='\0'; }
		informe<<key<<"="<<par<<endl;
	}
	fflush(stderr);
}
/*****************************************************************************
Lee una string
 *****************************************************************************/
void ReadParStr(string s,const char key[20], string par)
{
	int l,r;
    string str=key;
	l = strlen(key);
	if (s==key)
	{
		par=s;
        cout<<str<<" = "<<par<<endl;
	}
}
/*****************************************************************************
Lee varios decimal
 *****************************************************************************/
void ReadParMultD(char *s,const char key[20],int num, int *par)
{

	int l,r,i;
	char* pch;
	l = strlen(key);
	if (!strncmp(s,key,l))
	{
//		cout<<"ReadParMultD1 "<<s<<"   "<<key<<"   "<<num<<endl;
		pch = strtok (s," ");
		if (pch==NULL) { fprintf(stderr,"ERROR: No puede leerse %s en el fichero de parametros\n",key); exit(1); }
		i=0;
		while (pch != NULL)
		{
			pch = strtok (NULL, " ");
			par[i]=atoi(pch);
//			cout<<"pch: "<<pch<<endl;
//			cout<<"par[i]: "<<par[i]<<endl;
			i++;
//			cout<<"i: "<<i<<endl;
			if(i>num) break;
		}
	}
	fflush(stderr);
}
/*****************************************************************************
Lee  varios float
 *****************************************************************************/
void ReadParMultF(char *s,const char key[20],int num, double *par)
{
	int l,r,i;
	char* pch;
	l = strlen(key);
	if (!strncmp(s,key,l))
	{
		pch = strtok (s," ");
		if (pch==NULL) { fprintf(stderr,"ERROR: No puede leerse %s en el fichero de parametros\n",key); exit(1); }
		i=0;
		while (pch != NULL)
		{
			pch = strtok (NULL, " ");
			par[i]=atof(pch);
			i++;
		}
	}
	fflush(stderr);
}

int LeePotencialesOpticos(char *s,const char key[100],potencial_optico* pot,FILE* fp)
{
	int l,r,i,l2,l3;
	const char fin[]="FinPotencialesOpticos";
	const char flag[]="***************";
	char aux[500]="\0";
	l = strlen(key);
	l2 = strlen(fin);
	l3 = strlen(flag);
	i=0;
	if (!strncmp(s,key,l))
	{
		while(strncmp(aux,fin,l2))
		{
			fgets(aux,500,fp);
			ReadParD(aux,"id",&(pot[i].id));
			ReadParF(aux,"RealVolumen",&(pot[i].V));
			ReadParF(aux,"ImaginarioVolumen",&(pot[i].W));
			ReadParF(aux,"RealSpinOrbita",&(pot[i].Vso));
			ReadParF(aux,"ImaginarioSpinOrbita",&(pot[i].Wso));
			ReadParF(aux,"ImaginarioSuperficie",&(pot[i].Wd));
			ReadParF(aux,"RadioRealVolumen",&(pot[i].r0V));
			ReadParF(aux,"RadioImaginarioVolumen",&(pot[i].r0W));
			ReadParF(aux,"RadioCoulomb",&(pot[i].r0C));
			ReadParF(aux,"DifusividadRealVolumen",&(pot[i].aV));
			ReadParF(aux,"DifusividadImaginarioVolumen",&(pot[i].aW));
			ReadParF(aux,"RadioSpinOrbita",&(pot[i].rso));
			ReadParF(aux,"DifusividadSpinOrbita",&(pot[i].aso));
			ReadParF(aux,"RadioImaginarioSuperficie",&(pot[i].rWd));
			ReadParF(aux,"DifusividadImaginarioSuperficie",&(pot[i].aWd));
			if(!strncmp(aux,flag,3)) i++;
		}
		cout<<"Leidos "<<i<<" potenciales opticos"<<endl;
		fflush(stdout);
		return i;
	}
	return 0;
}

/*****************************************************************************
Lee los potenciales de campo medio del fichero de par�metros
 *****************************************************************************/
int LeePotencialesCampoMedio(char *s,const char key[100],potencial* pot,FILE* fp)
{
	int l,r,i,l2,l3;
	const char fin[]="FinCampoMedio";
	const char flag[]="***************";
	char aux[500]="\0";
	l = strlen(key);
	l2 = strlen(fin);
	i=0;
	if (!strncmp(s,key,l))
	{
		while(strncmp(aux,fin,l2))
		{
			fgets(aux,500,fp);
			ReadParD(aux,"id",&(pot[i].id));
			ReadParS(aux,"tipo",pot[i].tipo);
			ReadParF(aux,"VV",&(pot[i].V));
			ReadParF(aux,"VSO",&(pot[i].VSO));
			ReadParF(aux,"RV",&(pot[i].RV));
			ReadParF(aux,"aV",&(pot[i].aV));
			ReadParF(aux,"RSO",&(pot[i].RSO));
			ReadParF(aux,"aSO",&(pot[i].aSO));
			ReadParF(aux,"k",&(pot[i].k));
			ReadParF(aux,"rhc",&(pot[i].rhc));
			ReadParS(aux,"potfile",pot[i].file);
			if(!strncmp(aux,flag,3)) i++;
		}
		cout<<"Leidos "<<i<<" potenciales de campo medio"<<endl;
		fflush(stdout);
		return i;
	}
	return 0;
}
/*****************************************************************************
Lee los potenciales de campo medio del fichero de par�metros
 *****************************************************************************/
int LeeEstados(char *s,const char key[100],estado* st,FILE* fp)
{
	int l,r,i,l2,l3;
	const char fin[]="FinEstados";
	const char flag[]="***************";
	char aux[500]="\0";
	l = strlen(key);
	l2 = strlen(fin);
	i=0;
	if (!strncmp(s,key,l))
	{
		while(strncmp(aux,fin,l2))
		{
			fgets(aux,500,fp);
			ReadParD(aux,"id",&(st[i].id));
			ReadParD(aux,"l",&(st[i].l));
			ReadParF(aux,"j",&(st[i].j));
			ReadParD(aux,"nodos",&(st[i].nodos));
			ReadParF(aux,"spec",&(st[i].spec));
			ReadParF(aux,"energia",&(st[i].energia));
			ReadParS(aux,"file",st[i].file);
			if(!strncmp(aux,flag,3)) i++;
		}
		cout<<"Leidos "<<i<<" estados monoparticulares"<<endl;
		fflush(stdout);
		return i;
	}
	return 0;
}
/*****************************************************************************
Read Green's function
 *****************************************************************************/
int ReadGF(const char *fname,complejo** GF,double* r,int dimension,
	   double* einitial,double efinal,int ll,int dj)
{
  char aux[500];
  char* flag;
  const char fin[]="Ene";
  int l2=3;
  int n,m,l,cont1,cont2,jint,points;
  FILE *fp;
  complejo val,val2;
  float ene,j,r1,r2,RealPart,ImaginaryPart,r1old;
  fp = fopen(fname,"r");
  if (!fp) Error("Error opening Green's function file in ReadGF");
  n=0;
  m=0;
  flag=fgets(aux,500,fp);
  while(flag!=NULL)
    {
      sscanf(aux,"%*s %g %*s %d %*s %d %*s %d",&ene,&l,&jint,&points);
      j=jint/2.;
      flag=fgets(aux,500,fp);
      sscanf(aux,"%g %g %g %g",&r1,&r2,&RealPart,&ImaginaryPart);
      GF[0][0]=double(RealPart)+I*double(ImaginaryPart);
      r[0]=r1;
      cont1=0;
      cont2=0;
      r1old=r[0];
      flag=fgets(aux,500,fp);	
      while(strncmp(aux,fin,l2))
	{
	  if(cont1==0) {
	    r[cont2]=r2;
	  }
	  sscanf(aux,"%g %g %g %g",&r1,&r2,&RealPart,&ImaginaryPart);
	  if(r1==r1old){
	    cont2++;
	  }
	  else{
	    cont1++;
	    cont2=0;
	    r1old=r1;
	  }
	  GF[cont1][cont2]=double(RealPart)+I*double(ImaginaryPart);
	  //if(cont2>dimension) Error("Green's function file too big in ReadGF");  
	  flag=fgets(aux,500,fp);
	}
      if((ene>*einitial)&&(ene<efinal)&&(ll=l))
	{
	  cout<<"Energy: "<<ene<<"    L: "<<l<<"    j: "<<j<<"    points: "<<points<<endl;
	  // misc2<<"Energy: "<<ene<<"    L: "<<l<<"    j: "<<j<<"    points: "<<points<<endl;
	  // for(n=0;n<200;n++)
	  //   {
	  //     for(m=0;m<200;m++)
	  // 	{
	  // 	  misc2<<r[n]<<"  "<<r[m]<<"  "<<real(GF[n][m])<<"  "<<imag(GF[n][m])<<endl;
	  // 	}
	  //   }
	  *einitial=ene;
	  return 1;
	}
    }
  return 0;
}
/*****************************************************************************
Read non-local potential
 *****************************************************************************/
int ReadNLpot(const char *fname,complejo** potential,double* r,int dimension,
	   double energy,int ll,int dj)
{
  char aux[500];
  char* flag;
  const char fin[]="Ene";
  int l2=3;
  int n,m,l,cont1,cont2,jint,points;
  FILE *fp;
  complejo val,val2;
  float ene,j,r1,r2,RealPart,ImaginaryPart,r1old;
  fp = fopen(fname,"r");
  if (!fp) Error("Error opening potential file file in ReadNLpot");
  n=0;
  m=0;
  flag=fgets(aux,500,fp);
  while(flag!=NULL)
    {
      sscanf(aux,"%*s %g %*s %d %*s %d %*s %d",&ene,&l,&jint,&points);
      j=jint/2.;
      flag=fgets(aux,500,fp);
      sscanf(aux,"%g %g %g %g",&r1,&r2,&RealPart,&ImaginaryPart);
      potential[0][0]=double(RealPart)+I*double(ImaginaryPart);
      r[0]=r1;
      cont1=0;
      cont2=0;
      r1old=r[0];
      flag=fgets(aux,500,fp);	
      while(strncmp(aux,fin,l2))
	{
	  if(cont1==0) {
	    r[cont2]=r2;
	  }
	  sscanf(aux,"%g %g %g %g",&r1,&r2,&RealPart,&ImaginaryPart);
	  if(r1==r1old){
	    cont2++;
	  }
	  else{
	    cont1++;
	    cont2=0;
	    r1old=r1;
	  }
	  potential[cont1][cont2]=double(RealPart)+I*double(ImaginaryPart);
	  //if(cont2>dimension) Error("potential file too big in ReadNLpot");  
	  flag=fgets(aux,500,fp);
	}
      if((ene==energy)&&(l==ll)&&(jint==dj))
	{
	  //cout<<"Energy: "<<ene<<"    L: "<<l<<"    j: "<<j<<"    points: "<<points<<endl;
	  // misc3<<"Energy: "<<ene<<"    L: "<<l<<"    j: "<<j<<"    points: "<<points<<endl;
	  // for(n=0;n<200;n++)
	  //   {
	  //     for(m=0;m<200;m++)
	  // 	{
	  // 	  misc3<<r[n]<<"  "<<r[m]<<"  "<<real(potential[n][m])<<"  "<<imag(potential[n][m])<<endl;
	  // 	}
	  //   }
	  return 1;
	}
    }
  return 0;
}

void KoningDelaroche(double E,double N,double Z,double r,complejo* potencial_p,complejo* potencial_n,
		int l,double j,potencial_optico* pot_p,potencial_optico* pot_n)
{
	double A=N+Z;
	double v1p,v2p,v3p,v4p,wp1,wp2,dp1,dp2,dp3,vpso1,vpso2,wpso1,wpso2,rC,Vc,
	Epf,Vv,Wv,rV,aV,Wd,rD,aD,rSO,aSO,rwd,awd,radioV,radioW,radioWd,radio_coul,
	derivada,Wdd,delta_E,Vso,Wso,radioSO,ls,alpha,DVc;
	ofstream fp("KDPotential.txt");
	fp<<"Lab Energy: "<<E<<" MeV. N: "<<N<<". Z: "<<Z<<endl;
	fp<<"*********************************************"<<endl<<endl;
	ls=(j*(j+1.)-l*(l+1.)-0.75);
	alpha=(N-Z)/A;
	v1p=59.30+(21.*(N-Z)/A)-0.024*A;
	v2p=0.007067+4.23e-6*A;
	v3p=1.729e-5+1.136e-8*A;
	v4p=7e-9;
	wp1=14.667+0.009629*A;
	wp2=73.55+0.0795*A;
	dp1=16.0+16.*(N-Z)/A;
	dp2=0.0180+0.003802/(1+ exp((A-156.)/8.));
	dp3=11.5;
	vpso1=5.922+0.0030*A;
	vpso2=0.0040;
	wpso1=-3.1;
	wpso2=160.;
	rC=1.198+0.697*pow(A,(-2./3.))+12.994*pow(A,(-5./3.));
	Vc=(1.73/rC)*Z*pow(A,-0.333333333333);
	Epf=-8.4075 + 0.01378*A;
	Vv=v1p*(1.-v2p*(E-Epf)+v3p*(E-Epf)*(E-Epf)-v4p*(E-Epf)*(E-Epf)*(E-Epf))
			+Vc*v1p*(v2p-2.*v3p*(E-Epf)+3.*v4p*(E-Epf)*(E-Epf));
	Wv=wp1*(E-Epf)*(E-Epf)/((E-Epf)*(E-Epf)+wp2*wp2);
	rV=1.3039-0.4054*pow(A,(-1./3.));
	aV=0.6778-1.487e-4*A;
	Wd=exp(-dp2*(E-Epf))*dp1*(E-Epf)*(E-Epf)/((E-Epf)*(E-Epf)+dp3*dp3);
	rD=1.3424-0.01585*pow(A,(1./3.));
	aD=0.5187+5.205e-4*A;
	rwd=1.3424-0.01585*pow(A,(1./3.));
	awd=0.5187+5.205e-4*A;
	Vso=vpso1*exp(-vpso2*(E-Epf));
	Wso=wpso1*(E-Epf)*(E-Epf)/((E-Epf)*(E-Epf)+wpso2*wpso2);
	radioSO=1.1854-0.647*pow(A,-1./3.);
	aSO=0.59;
	radioV=rV*pow(A,0.33333333333333);
	radioSO=radioSO*pow(A,0.33333333333333);
	radioW=rV*pow(A,0.33333333333333);
	radioWd=rwd*pow(A,0.33333333333333);
	radio_coul=rC*pow(A,0.33333333333333);
//	cout<<"Vv: "<<Vv<<endl;
	if(j==0.)
		*potencial_p=-Vv/(1.+exp((r-radioV)/aV))-I*Wv/(1.+exp((r-radioW)/aV))
		-4.*I*Wd*exp((r-radioWd)/awd)/((1.+exp((r-radioWd)/awd))*(1.+exp((r-radioWd)/awd)));
	else
		*potencial_p=-Vv/(1.+exp((r-radioV)/aV))-I*Wv/(1.+exp((r-radioW)/aV))
		-4.*I*Wd*exp((r-radioWd)/awd)/((1.+exp((r-radioWd)/awd))*(1.+exp((r-radioWd)/awd)))
		-2.*(ls*Vso)*exp((r-radioSO)/aSO)/(r*aSO*(1.+exp((r-radioSO)/aSO))*(1.+exp((r-radioSO)/aSO)));
//   misc2<<r<<"  "<<real(*potencial_p)<<"  "<<imag(*potencial_p)<<endl;
//	cout<<" rV: "<<rV<<" aV: "<<aV<<" v1p: "<<v1p<<" v2p: "<<v2p<<" v3p: "<<v3p<<" Epf: "<<Epf<<" Vv: "<<Vv<<endl;

	pot_p->V=Vv;
	pot_p->Vso=2.*Vso;
	pot_p->W=Wv;
	pot_p->Wd=Wd;
	pot_p->aV=aV;
	pot_p->aW=aV;
	pot_p->aWd=awd;
	pot_p->aso=aSO;
	pot_p->r0C=rC;
	pot_p->r0V=rV;
	pot_p->r0W=rV;
	pot_p->rWd=rwd;
	pot_p->radioV=radioV;
	pot_p->radioso=radioSO;
	pot_p->radioWd=radioWd;
	pot_p->rso=radioWd;
	pot_p->radio_coul=radio_coul;

	fp<<"Protons: "<<endl<<endl;
	fp<<"RealVolumen "<<Vv<<endl<<
			"ImaginarioVolumen  "<<Wv<<endl<<
			"RealSpinOrbita  "   <<Vso<<endl<<
			"ImaginarioSpinOrbita	"<<0.<<endl<<
			"ImaginarioSuperficie  " <<Wd<<endl<<
			"RadioRealVolumen  "   <<rV<<endl<<
			"RadioCoulomb  "            <<rC<<endl<<
			"RadioImaginarioVolumen  "       <<rV<<endl<<
			"DifusividadRealVolumen  "         <<aV<<endl<<
			"DifusividadImaginarioVolumen  "   <<aV<<endl<<
			"RadioSpinOrbita    "         	<<radioSO<<endl<<
			"DifusividadSpinOrbita  "       <<aSO<<endl<<
			"RadioImaginarioSuperficie  "          <<rwd<<endl<<
			"DifusividadImaginarioSuperficie "    <<awd<<endl;


	v1p=59.30-21.*(N-Z)/A-0.024*A;
	v2p=0.007228-1.48e-6*A;
	v3p=1.994e-5-2.e-8*A;
	v4p=7.e-9;
	wp1=12.195+0.0167*A;
	wp2=73.55+0.0795*A;
	dp1=16.0-16*(N-Z)/A;
	dp2=0.0180+0.003802/(1.+ exp((A-156.)/8.));
	dp3=11.5;
	vpso1=5.922+0.0030*A;
	vpso2=0.0040;
	wpso1=-3.1;
	wpso2=160.;
	rC=1.198+0.697*pow(A,(-2./3.))+12.994*pow(A,(-5./3.));
	Vc=1.73/(rC*Z*pow(A,-0.333333333333));
	Epf=-11.2814+0.02646*A;

	Vv=v1p*(1.-v2p*(E-Epf)+v3p*(E-Epf)*(E-Epf)-v4p*(E-Epf)*(E-Epf)*(E-Epf));
	Wv=wp1*(E-Epf)*(E-Epf)/((E-Epf)*(E-Epf)+wp2*wp2);
	rV=1.3039-0.4054*pow(A,(-1./3.));
	aV=0.6778-1.487e-4*A;
	Wd=exp(-dp2*(E-Epf))*dp1*(E-Epf)*(E-Epf)/((E-Epf)*(E-Epf)+dp3*dp3);
	rD=1.3424-0.01585*pow(A,(1./3.));
	aD=0.5446-1.656e-4*A;
	rwd=1.3424-0.01585*pow(A,(1./3.));
	awd=0.5446-1.656e-4*A;
	Vso=vpso1*exp(-vpso2*(E-Epf));
	Wso=wpso1*(E-Epf)*(E-Epf)/((E-Epf)*(E-Epf)+wpso2*wpso2);
	radioSO=1.1854-0.647*pow(A,-1./3.);
	aSO=0.59;




	radioV=rV*pow(A,0.33333333333333);
	radioSO=radioSO*pow(A,0.33333333333333);
	radioW=rV*pow(A,0.33333333333333);
	radioWd=rwd*pow(A,0.33333333333333);
	radio_coul=rC*pow(A,0.33333333333333);


    /////////////////////////////////////////////////
	// For 95Mo:
//	Vv=50.;
//	if(Wd<=4.) Wd=4.;
	////////////////////////////////////////////////

    /////////////////////////////////////////////////
	// For 93Nb:
//	Vv=50.3;
//	if(Wd<=4.) Wd=4.;
	////////////////////////////////////////////////

//	cout<<" rV: "<<rV<<" aV: "<<aV<<" v1p: "<<v1p<<" v2p: "<<v2p<<" v3p: "<<v3p<<" Epf: "<<Epf<<endl;
	if(j==0.)
	*potencial_n=-Vv/(1.+exp((r-radioV)/aV))-I*Wv/(1.+exp((r-radioW)/aV))
			-4.*I*Wd*exp((r-radioWd)/awd)/((1.+exp((r-radioWd)/awd))*(1.+exp((r-radioWd)/awd)));
	else
	*potencial_n=-Vv/(1.+exp((r-radioV)/aV))-I*Wv/(1.+exp((r-radioW)/aV))
			-4.*I*Wd*exp((r-radioWd)/awd)/((1.+exp((r-radioWd)/awd))*(1.+exp((r-radioWd)/awd)))
			-2.*(ls*Vso)*exp((r-radioSO)/aSO)/(r*aSO*(1.+exp((r-radioSO)/aSO))*(1.+exp((r-radioSO)/aSO)));

//	misc3<<r<<"  "<<real(*potencial_n)<<"  "<<imag(*potencial_n)<<endl;
//	cout<<"Potencial Koning-Delaroche. Radio reducido: "<<rV<<",   radio: "<<radioV<<",   difusividad: "<<aV<<endl;
//	cout<<"Radio spin-orbita: "<<radioSO<<",   profundidad spin orbita: "<<Vso<<",   difusividad spin orbita: "<<aSO<<endl;

	pot_n->V=Vv;
	pot_n->Vso=Vso;
	pot_n->W=Wv;
	pot_n->Wd=Wd;
	pot_n->aV=aV;
	pot_n->aW=aV;
	pot_n->aWd=awd;
	pot_n->aso=aSO;
	pot_n->r0C=rC;
	pot_n->r0V=rV;
	pot_n->r0W=rV;
	pot_n->rWd=rwd;
	pot_n->radioV=radioV;
	pot_n->radioso=radioSO;
	pot_n->radioWd=radioWd;
	pot_n->rso=radioWd;

	fp<<"Neutrons: "<<endl<<endl;
	fp<<"RealVolumen "<<Vv<<endl<<
			"ImaginarioVolumen  "<<Wv<<endl<<
			"RealSpinOrbita  "   <<Vso<<endl<<
			"ImaginarioSpinOrbita	"<<0.<<endl<<
			"ImaginarioSuperficie  " <<Wd<<endl<<
			"RadioRealVolumen  "   <<rV<<endl<<
			"RadioCoulomb  "            <<rC<<endl<<
			"RadioImaginarioVolumen  "       <<rV<<endl<<
			"DifusividadRealVolumen  "         <<aV<<endl<<
			"DifusividadImaginarioVolumen  "   <<aV<<endl<<
			"RadioSpinOrbita    "         	<<radioSO<<endl<<
			"DifusividadSpinOrbita  "       <<aSO<<endl<<
			"RadioImaginarioSuperficie  "          <<rwd<<endl<<
			"DifusividadImaginarioSuperficie "    <<awd<<endl;
//	misc1<<E<<"  "<<Vv<<"  "<<Wd<<"  "<<Wv<<"  "<<Vso<<"  "<<rV<<"  "<<aV<<"  "<<rwd<<"  "<<awd<<endl;
//cout<<" radio: "<<Vv/(1.+exp((r-radioV)/aV))<<"   "<<Vv<<"   "<<r<<"   "<<radioV<<"   "<<1.+exp((r-radioV)/aV)<<endl;
//	*potencial_n=-2.*(ls*Vso)*exp((r-radioSO)/aSO)/(r*aSO*(1.+exp((r-radioSO)/aSO))*(1.+exp((r-radioSO)/aSO)));

	/////////////////////////////// Calculo derivada absorcion de superficie  //////////////////////////
//	delta_E=0.001;
//	E+=delta_E;
//	Wdd=v1p*(1-v2p*(E-Epf)+v3p*(E-Epf)*(E-Epf)-v4p*(E-Epf)*(E-Epf)*(E-Epf));
//	derivada=(Wdd-Vv)/delta_E;
//	misc1<<E-delta_E<<"  "<<Wd<<"  "<<derivada<<"  "<<Wdd<<"  "<<Vv<<endl;
}
void CH89(double E,double N,double Z,double r,complejo* potencial_p,complejo* potencial_n,
		int l,double j,potencial_optico* pot_p,potencial_optico* pot_n)
{
	double A=N+Z;
	double V0,Vt,Ve,r0,r00,a0,rc,rc0,Vso,rso,rso0,aso,Wv0,Wve0,Wvew,Ws0,
	      Wst,Wse0,Wsew,rW,rW0,aW,Vrp,Vrn,R0,Rc,Rso,Wvn,Wvp,Wsp,Wsn,Rw,Ecp,Ecn,ls;
	ofstream fp("CH89Potential.txt");
	fp<<"Lab Energy: "<<E<<" MeV. N: "<<N<<". Z: "<<Z<<endl;
	fp<<"*********************************************"<<endl<<endl;
	ls=(j*(j+1.)-l*(l+1.)-0.75);
	V0=52.9;
	Vt=13.1;
	Ve=-0.299;
	r0=1.25;
	r00=-0.225;
	a0=0.69;
	rc=1.238;
	rc0=0.116;
	Vso=5.9;
	rso=1.34;
	rso0=-1.2;
	aso=0.63;
	Wv0=7.8;
	Wve0=35.;
	Wvew=16.;
	Ws0=10.;
    Wst=18.;
    Wse0=36.;
    Wsew=37.;
    rW=1.33;
    rW0=-0.42;
    aW=0.69;
    Rc=rc*pow(A,0.3333333333)+rc0;
    Ecp=1.73*Z/Rc;
    Ecn=0;
    Vrp=V0+(Vt*(N-Z)/A)+(E-Ecp)*Ve;
    Vrn=V0-(Vt*(N-Z)/A)+(E-Ecn)*Ve;
    R0=r0*pow(A,0.3333333333)+r00;
    Rso=rso*pow(A,0.3333333333)+rso0;
    Wvn=Wv0/(1.+exp((Wve0-(E-Ecn))/(Wvew)));
    Wvp=Wv0/(1.+exp((Wve0-(E-Ecp))/(Wvew)));
    Wsp=(Ws0+Wst*((N-Z)/A))/(1.+exp(((E-Ecp)-Wse0)/(Wsew)));
    Wsn=(Ws0-Wst*((N-Z)/A))/(1.+exp(((E-Ecn)-Wse0)/(Wsew)));
    Rw=rW*pow(A,0.3333333333)+rW0;


	pot_n->V=Vrn;
	pot_n->Vso=Vso;
	pot_n->W=Wvn;
	pot_n->Wd=Wsn;
	pot_n->aV=a0;
	pot_n->aW=aW;
	pot_n->aWd=aW;
	pot_n->aso=aso;
	pot_n->r0C=Rc;
	pot_n->r0V=R0;
	pot_n->r0W=Rw;
	pot_n->rWd=Rw;
	pot_n->radioV=R0;
	pot_n->radioso=Rso;
	pot_n->radioWd=Rw;
	pot_n->rso=rso;

    /////////////////////////////////////////////////
	// For 95Mo:
//	Vrn=50.;
//	if(Wvn<=4.) Wvn=4.;
	////////////////////////////////////////////////

	if(j==0.)
		*potencial_n=-Vrn/(1.+exp((r-R0)/a0))-I*Wvn/(1.+exp((r-Rw)/aW))
		-4.*I*Wsn*exp((r-Rw)/aW)/((1.+exp((r-Rw)/aW))*(1.+exp((r-Rw)/aW)));
	else
		*potencial_n=-Vrn/(1.+exp((r-R0)/a0))-I*Wvn/(1.+exp((r-Rw)/aW))
		-4.*I*Wsn*exp((r-Rw)/aW)/((1.+exp((r-Rw)/aW))*(1.+exp((r-Rw)/aW)))
		-2.*(ls*Vso)*exp((r-Rso)/aso)/(r*aso*(1.+exp((r-Rso)/aso))*(1.+exp((r-Rso)/aso)));

	fp<<"Neutrons: "<<endl<<endl;
	fp<<"RealVolumen "<<Vrn<<endl<<
			"ImaginarioVolumen  "<<Wvn<<endl<<
			"RealSpinOrbita  "   <<Vso<<endl<<
			"ImaginarioSpinOrbita	"<<0.<<endl<<
			"ImaginarioSuperficie  " <<Wsn<<endl<<
			"RadioRealVolumen  "   <<r0<<endl<<
			"RadioCoulomb  "            <<rc<<endl<<
			"RadioImaginarioVolumen  "       <<rW<<endl<<
			"DifusividadRealVolumen  "         <<a0<<endl<<
			"DifusividadImaginarioVolumen  "   <<aW<<endl<<
			"RadioSpinOrbita    "         	<<rso<<endl<<
			"DifusividadSpinOrbita  "       <<aso<<endl<<
			"RadioImaginarioSuperficie  "          <<rW<<endl<<
			"DifusividadImaginarioSuperficie "    <<aW<<endl;


	pot_p->V=Vrp;
	pot_p->Vso=Vso;
	pot_p->W=Wvp;
	pot_p->Wd=Wsp;
	pot_p->aV=a0;
	pot_p->aW=aW;
	pot_p->aWd=aW;
	pot_p->aso=aso;
	pot_p->r0C=Rc;
	pot_p->r0V=R0;
	pot_p->r0W=Rw;
	pot_p->rWd=Rw;
	pot_p->radioV=R0;
	pot_p->radioso=Rso;
	pot_p->radioWd=Rw;
	pot_p->rso=rso;



	if(j==0.)
		*potencial_p=-Vrp/(1.+exp((r-R0)/a0))-I*Wvp/(1.+exp((r-Rw)/aW))
		-4.*I*Wvp*exp((r-Rw)/aW)/((1.+exp((r-Rw)/aW))*(1.+exp((r-Rw)/aW)));
	else
		*potencial_n=-Vrn/(1.+exp((r-R0)/a0))-I*Wvp/(1.+exp((r-Rw)/aW))
		-4.*I*Wsp*exp((r-Rw)/aW)/((1.+exp((r-Rw)/aW))*(1.+exp((r-Rw)/aW)))
		-2.*(ls*Vso)*exp((r-Rso)/aso)/(r*aso*(1.+exp((r-Rso)/aso))*(1.+exp((r-Rso)/aso)));
	fp<<"Protons: "<<endl<<endl;
	fp<<"RealVolumen "<<Vrp<<endl<<
			"ImaginarioVolumen  "<<Wvp<<endl<<
			"RealSpinOrbita  "   <<Vso<<endl<<
			"ImaginarioSpinOrbita	"<<0.<<endl<<
			"ImaginarioSuperficie  " <<Wsp <<endl<<
			"RadioRealVolumen  "   <<r0<<endl<<
			"RadioCoulomb  "            <<rc<<endl<<
			"RadioImaginarioVolumen  "       <<rW<<endl<<
			"DifusividadRealVolumen  "         <<a0<<endl<<
			"DifusividadImaginarioVolumen  "   <<aW<<endl<<
			"RadioSpinOrbita    "         	<<rso<<endl<<
			"DifusividadSpinOrbita  "       <<aso<<endl<<
			"RadioImaginarioSuperficie  "          <<rW<<endl<<
			"DifusividadImaginarioSuperficie "    <<aW<<endl;
}
///////////////////////////////////////////////////////////////////////
//                                                                   //
//       Potencial optico de Han, Shi and Shen                        //
//                                                                   //
//////////////////////////////////////////////////////////////////////
void HanShiShen(double E,double N,double Z)
{
	double A=N+Z;
	double V0,V1,V2,V3,V4,W0,W1,W2,U0,U1,U2,U3;
	double VR,WD,WS,WSO,VSO,RR,aR,RD,aD,RS,aS,RSO,aSO,RC,ad,as;
	double rR,rD,rS,rSO,rC;
	ofstream fp("HanShiShenPotential.txt");
	fp<<"Lab Energy: "<<E<<" MeV. N: "<<N<<". Z: "<<Z<<endl;
	fp<<"*********************************************"<<endl<<endl;
	V0=82.18;
	V1=-0.148;
	V2=-0.000886;
	V3=-34.811;
	V4=1.058;
	W0=20.968;
	W1=-0.0794;
	W2=-43.398;
	U0=-4.916;
	U1=0.0555;
	U2=0.0000442;
	U3=35.0;
	VSO=3.703;
	WSO=-0.206;
	rR=1.174;
	rD=1.328;
	rS=1.563;
	rSO=1.234;
	rC=1.698;
	aR=0.809;
	ad=0.465;
	as=0.7;
	aD=ad+0.045*pow(A,1./3.);
	aS=as+0.045*pow(A,1./3.);
	aSO=0.813;
	VR=V0+V1*E+V2*E*E+V3*((N-Z)/A)+V4*Z/pow(A,1./3.);
	WD=W0+W1*E+W2*(N-Z)/A;
	WS=U0+U1*E+U2*E*E+U3*((N-Z)/A);
	if(WS<0.) WS=0.;
	fp<<"RealVolumen "<<VR<<endl<<
			"ImaginarioVolumen  "<<WS<<endl<<
			"RealSpinOrbita  "   <<VSO<<endl<<
			"ImaginarioSpinOrbita	"<<WSO<<endl<<
			"ImaginarioSuperficie  " <<WD<<endl<<
			"RadioRealVolumen  "   <<rR<<endl<<
			"RadioCoulomb  "            <<rC<<endl<<
			"RadioImaginarioVolumen  "       <<rS<<endl<<
			"DifusividadRealVolumen  "         <<aR<<endl<<
			"DifusividadImaginarioVolumen  "   <<aS<<endl<<
			"RadioSpinOrbita    "         	<<rSO<<endl<<
			"DifusividadSpinOrbita  "       <<aSO<<endl<<
			"RadioImaginarioSuperficie  "          <<rD<<endl<<
			"DifusividadImaginarioSuperficie "    <<aD<<endl;
}
void GeneraPotencialOptico(struct parametros *parm,struct potencial_optico *potencial,double m1,double m2)
{
	int n;
	double delta_r;
	delta_r=parm->radio/double(parm->puntos);
	if(m1<1. || m2<1.) Error("Masa menor de 1");
	if(m1>3. && m2>3.)
	{
		potencial->radioV=potencial->r0V*(pow(m1,0.33333333333333)+pow(m2,0.33333333333333));
		potencial->radioW=potencial->r0W*(pow(m1,0.33333333333333)+pow(m2,0.33333333333333));
		potencial->radioso=potencial->rso*(pow(m1,0.33333333333333)+pow(m2,0.33333333333333));
		potencial->radioWd=potencial->rWd*(pow(m1,0.33333333333333)+pow(m2,0.33333333333333));
		potencial->radio_coul=potencial->r0C*(pow(m1,0.33333333333333)+pow(m2,0.33333333333333));
	}
	if(m1>3. && m2<=3.)
	{
		potencial->radioV=potencial->r0V*pow(m1,0.33333333333333);
		potencial->radioW=potencial->r0W*pow(m1,0.33333333333333);
		potencial->radioso=potencial->rso*pow(m1,0.33333333333333);
		potencial->radioWd=potencial->rWd*pow(m1,0.33333333333333);
		potencial->radio_coul=potencial->r0C*pow(m1,0.33333333333333);
	}
	if(m1<=3. && m2>3.)
	{
		potencial->radioV=potencial->r0V*pow(m2,0.33333333333333);
		potencial->radioW=potencial->r0W*pow(m2,0.33333333333333);
		potencial->radioso=potencial->rso*pow(m2,0.33333333333333);
		potencial->radioWd=potencial->rWd*pow(m2,0.33333333333333);
		potencial->radio_coul=potencial->r0C*pow(m2,0.33333333333333);
	}
	for(n=0;n<parm->puntos;n++)
	{
		potencial->r[n]=delta_r*(n+1.);
		potencial->pot[n]=-potencial->V/(1.+exp((potencial->r[n]-potencial->radioV)/potencial->aV))-I*potencial->W/
				(1.+exp((potencial->r[n]-potencial->radioW)/potencial->aW))-4.*I*potencial->Wd*
				exp((potencial->r[n]-potencial->radioWd)/potencial->aWd)/((1.+exp((potencial->r[n]-potencial->radioWd)/potencial->aWd))
						*(1.+exp((potencial->r[n]-potencial->radioWd)/potencial->aWd)));
	}
	potencial->puntos=parm->puntos;
}
void GeneraPotencialCM(struct parametros *parm,struct potencial *potencial)
{
	int i;
	double delta_r;
	delta_r=parm->radio/double(parm->puntos);
	potencial->puntos=parm->puntos;
	if(!strcmp(potencial->tipo,"ws"))
	{
		for (i=0;i<potencial->puntos;i++) {
			potencial->r[i]=delta_r*(i+1);
			potencial->pot[i]=-(potencial->V)/(1.+exp((potencial->r[i]-potencial->RV)/potencial->aV));
		}
	}
	if(!strcmp(potencial->tipo,"tang"))
	{
		for (i=0;i<potencial->puntos;i++) {
			potencial->r[i]=delta_r*(i+1);
			if(potencial->r[i]>potencial->rhc) potencial->pot[i]=-potencial->V*exp(-potencial->k*(potencial->r[i]-potencial->rhc));
			if (potencial->r[i]<=potencial->rhc) potencial->pot[i]=potencial->V;
		}
	}
}
/***************************************************************
 * Genera coeficientes de Clebsh-Gordan <l1 m1 l2 m2|J M>
 ***************************************************************/
double ClebsGordan(float l1,float m1,float l2,float m2,float J,float M)
{
	double cg;
	if ((abs(l1-l2)>J) || (l1+l2<J)) return 0.;
	if (abs(m1)>l1) return 0.;
	if (abs(m2)>l2) return 0.;
	if (m1+m2!=M) return 0.;
	cg=pow(-1.,l1-l2+M)*sqrt(2.*J+1.)*gsl_sf_coupling_3j (int(2.*l1),int(2.*l2),int(2.*J),int(2.*m1),int(2.*m2),-int(2.*M));
	return cg;
}
/*****************************************************************************
Escribe las energias y las funciones de onda
 *****************************************************************************/
void EscribeEstados(int puntos,estado* st,int numero_estados,struct parametros *parm)
{
	ofstream fpen(parm->fl_energias);
	ofstream fpst(parm->fl_funondas);
	int n,i;
	fpen<<"     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<"\n";
	fpen<<"     +                                                      +"<<"\n";
	fpen<<"     +         Estados de part�cula independiente           +"<<"\n";
	fpen<<"     +                                                      +"<<"\n";
	fpen<<"     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<"\n"<<"\n"<<"\n";
	fpen<<"Indice   "<<"Nodos"<<"          L"<<"             J"<<"             Energia"<<"\n";
	for(n=0;n<numero_estados;n++){
		fpen<<"  "<<n<<"........"<<st[n].nodos<<"..........."<<" "<<st[n].l<<"..........."<<" "<<st[n].j<<"..........."<<st[n].energia<<"\n";
		cout<<"  "<<n<<"........"<<st[n].nodos<<"..........."<<" "<<st[n].l<<"..........."<<" "<<st[n].j<<"..........."<<st[n].energia<<"\n";
	}
	for(i=0;i<puntos;i++){
		fpst<<st[0].r[i]<<"  ";
		for(n=0;n<numero_estados;n++){
//			fpst<<real(st[n].wf[i])<<"  ";
			fpst<<real(st[n].wf[i]*st[0].r[i])<<"  "<<real(st[n].wf[i])<<"  ";
		}
		fpst<<endl;
	}
	fpen.close();
	fpst.close();
}
/*****************************************************************************
Escribe el potencial
 *****************************************************************************/
void EscribePotencial(int puntos,potencial* pot,int numero_potenciales,struct parametros *parm)
{
	ofstream fp(parm->fl_potcm);
	int n,i;
	cout<<"potenciales... "<<numero_potenciales<<endl;
	for(i=0;i<puntos;i++){
		fp<<pot[0].r[i]<<"  ";
		for(n=0;n<numero_potenciales;n++){
			fp<<pot[n].pot[i]<<"  ";
		}
		fp<<endl;
	}
	fp.close();
}
//////////////////////////////////////////////
//    Inicializador de matrices  complejo     //
//////////////////////////////////////////////
complejo** matriz_cmpx(int dim1,int dim2)
{
  int n,m;
  complejo** mat;
  mat=new complejo*[dim1];
  for(n=0;n<dim1;n++)
  {
    mat[n]=new complejo[dim2];
  }
  for(n=0;n<dim1;n++)
  {
    for(m=0;m<dim2;m++)
    {
      mat[n][m]=0.;
    }
  }
  return (mat);
}
//////////////////////////////////////////////
//    Inicializador de matrices  double     //
//////////////////////////////////////////////
double** matriz_dbl(int dim1,int dim2)
{
  int n,m;
  double** mat;
  mat=new double*[dim1];
  for(n=0;n<dim1;n++)
  {
    mat[n]=new double[dim2];
  }
  for(n=0;n<dim1;n++)
  {
    for(m=0;m<dim2;m++)
    {
      mat[n][m]=0.;
    }
  }
  return (mat);
}
/*****************************************************************************
Genera "regla" puntos de integracion gaussiana. "abs" contiene las posiciones
y "w" los pesos
 *****************************************************************************/
void GaussLegendre(double* absi,double* w,int regla)
{
  double inf=0.;
  char c[100],d[200];
  int n;
//  if(regla<=60)
//  {
//	  FILE *gl;
//	  sprintf(c,"%d",regla);
//	  strcpy(d,"C:\\Gregory\\workspace\\coops\\input\\GaussLegendre\\GL");
//	  strcat(d,c);
//	  strcat(d,".txt");
//	  gl=fopen(d,"r");
//	  for(n=0;n<regla;n++)
//	  {
//		  fscanf(gl,"%le",&inf);
//		  absi[n]=inf;
//		  fscanf(gl,"%le\n",&inf);
//		  w[n]=inf;
//	  }
//	  fclose(gl);
//  }
//  else
  {
	  LegendreRoots(regla,absi,w);
  }
}
void LegendreRoots(int regla, double* absi, double* w)
{
	int i;
	double x, x1;
	for (i=1;i<=regla;i++) {
		x=cos(PI*(i-.25)/(regla+.5));
		do {
			x1=x;
			x-=gsl_sf_legendre_Pl(regla,x)/ lege_diff(regla, x);
//			cout<<x<<"  "<<x1<<endl;
		} while (abs(x-x1)>1.e-10);
		absi[i-1]=-x;
		x1=lege_diff(regla,x);
		w[i-1]=2/((1-x*x)*x1*x1);
	}
}
double lege_diff(int n, double x)
{
//	cout<<"lege_diff1"<<endl;
	return n*(x *gsl_sf_legendre_Pl(n,x)-gsl_sf_legendre_Pl(n-1,x))/(x*x-1);
}
double NormalizaD(distorted_wave* st1,distorted_wave* st2, double radio, int pts, char s) {
	int regla_r, nr, indice;
	double ar, br, norma, r,radio_medio;
	double step=radio/double(pts);
	regla_r = 60;
	double* wr = new double[regla_r];
	double* absr = new double[regla_r];
	complejo sum = 0.;
	ar = 0.;
	br = radio;
	GaussLegendre(absr, wr, regla_r);
	for (nr = 0; nr < regla_r; nr++) {
		r = ar + (br - ar) * (absr[nr] + 1.) / 2.;
		indice = int(ceil(r / step)) - 1;
		if (indice > pts - 1) indice = pts - 1;
		sum+=st1->wf[indice]*st2->wf[indice]*r*r*wr[nr];
	}
	norma = abs(sum) * (br - ar) / 2.;
	radio_medio=0.;
	if (s == 's') {
		for (nr = 0; nr < pts; nr++) {
			r = (nr + 1) * step;
			st1->wf[nr] = st1->wf[nr] / sqrt(norma);
			radio_medio+=abs(st1->wf[nr]*st1->wf[nr])*r*r*r*step;
		}
	}
//	cout<<"Radio: "<<radio_medio<<endl;
	delete[] wr;
	delete[] absr;
	return norma;
}
/*****************************************************************************
Calcula el solapamiento entre "st1" y "st2". Si "s"=s, divide "st1" por la raiz
cuadrada de este valor.
 *****************************************************************************/
double Normaliza(estado* st1,estado* st2, double radio, int pts, char s) {
	int regla_r, nr, indice;
	double ar, br, norma, r,radio_medio;
	double step=radio/double(pts);
	regla_r = 60;
	double* wr = new double[regla_r];
	double* absr = new double[regla_r];
	complejo sum = 0.;
	ar = 0.;
	br = radio;
//	cout<<"Normaliza1"<<endl;
	GaussLegendre(absr, wr, regla_r);
	for (nr = 0; nr < regla_r; nr++) {
		r = ar + (br - ar) * (absr[nr] + 1.) / 2.;
		indice = int(ceil(r / step)) - 1;
		if (indice > pts - 1) indice = pts - 1;
		sum+=st1->wf[indice]*st2->wf[indice]*r*r*wr[nr];
	}
	norma = abs(sum) * (br - ar) / 2.;
	radio_medio=0.;
	if (s == 's') {
		for (nr = 0; nr < pts; nr++) {
			r = (nr + 1) * step;
			st1->wf[nr] = st1->wf[nr] / sqrt(norma);
			radio_medio+=abs(st1->wf[nr]*st1->wf[nr])*r*r*r*r*step;
		}
	}
	delete[] wr;
	delete[] absr;
	return norma;
}
/*****************************************************************************
Interpolacion lineal de una funcion compleja
 *****************************************************************************/
complejo interpola_cmpx(complejo* funcion,double* r,double posicion,int puntos)
{
	double delta_r;
	int idx;
	if (puntos<3) Error("El numero de puntos debe ser mayor que 3!");
	delta_r=r[puntos-1]-r[puntos-2];
	idx=int(ceil(posicion/delta_r));
//	misc1<<delta_r<<"  "<<idx<<"  "<<r[puntos-1]<<"  "<<r[puntos-2]<<endl;
	if(idx>puntos-2) return funcion[puntos-1];
	if(idx<1) return funcion[0];
	return (funcion[idx]+(funcion[idx+1]-funcion[idx])*(posicion-r[idx])/delta_r);
}
/*****************************************************************************
Interpolacion lineal de una funcion real
 *****************************************************************************/
double interpola_dbl(double* funcion,double* r,double posicion,int puntos)
{
	int indice;
	double a, b, f0, f1, f2, x0, x1, x2,delta_r;
	if (puntos<3) Error("El numero de puntos debe ser mayor que 3!");
	delta_r=r[puntos-1]-r[puntos-2];
	indice = int(ceil(posicion/delta_r)) - 1;
	if (indice > puntos - 2)
		return funcion[puntos-1];
	if (indice <= 0)
		return funcion[0];
	f0 = funcion[indice-1];
	f1 = funcion[indice];
	f2 = funcion[indice+1];
	x0 = delta_r * (indice);
	x1 = delta_r * (indice + 1.);
	x2 = delta_r * (indice + 2.);
	b = ((x2 - x0) * (x2 - x0) * (f1 - f0) + (x1 - x0) * (x1 - x0) * (f0 - f2))
			/ ((x2 - x0) * (x2 - x0) * (x1 - x0) + (x1 - x0) * (x1 - x0) * (x0
					- x2));
	a = (f2 - f0 - b * (x2 - x0)) / ((x2 - x0) * (x2 - x0));
	if (f0 == 0.)
		return 0.;
	return a * (posicion - x0) * (posicion - x0) + b * (posicion - x0) + f0;
}
complejo GeneraDWspin(distorted_wave* funcion,potencial_optico *v, double q1q2, double masa,double radio_max,
		int puntos,double radio_match,ofstream* fp)
{
	int i, i_1, i_2,status1,status2;
	double hbarx, dd, radio_1, radio_2, q, x, y, ex1, ex2, etac,delta_r,spinorbit;
	complejo delta, derivada_log, fu1, fu2, factor;
	complejo *potencial=new complejo[puntos];
	gsl_complex deltagsl;
	gsl_sf_result F1, G1, F2, G2, Fp, Gp;
	if(funcion[0].energia<=0.) Error("Energia negativa o 0 en GeneraDW");
	delta_r=radio_max/double(puntos);
	if(radio_match>radio_max-4.*delta_r) Error("Radio de matching demasiado grande en GeneraDW");
	hbarx=HC*HC/(2.*AMU*masa);
	dd=delta_r*delta_r/hbarx;
	q=sqrt(funcion[0].energia/hbarx);
	etac=q1q2*masa*E2HC*AMU/(HC*q);
	funcion->puntos=puntos;
	funcion->radio=radio_max;
//	cout<<"pot down"<<v_down->pot[3]<<endl;
	spinorbit =(funcion->j)*((funcion->j)+1.)-funcion->l*(funcion->l+1.)-(funcion->spin)*((funcion->spin)+1.); //T�rmino de spin-�rbita
	/* actualizacion del potencial con los t�rminos de Coulomb, centr�fugo y de spin-�rbita*/
	for (i=0;i<puntos-1;i++) {
		if(v->r[i]>=v->radio_coul) potencial[i]=v->pot[i]+E_CUADRADO*q1q2/v->r[i]+
				(funcion[0].l*(funcion[0].l+1.))*hbarx /(v->r[i]*v->r[i])
				-2.*spinorbit*v->Vso*exp((v->r[i]-v->radioso)/v->aso)
		/((v->aso*v->r[i])*(1.+exp((v->r[i]-v->radioso)/v->aso))*(1.+exp((v->r[i]-v->radioso)/v->aso)));
		if(v->r[i]<v->radio_coul) potencial[i]=v->pot[i]+E_CUADRADO*q1q2*(3.-(v->r[i]/v->radio_coul)* (v->r[i]/v->radio_coul))/(2.*v->radio_coul)+
				(funcion[0].l*(funcion[0].l+1.))*hbarx /(v->r[i]*v->r[i])
				-2.*spinorbit*v->Vso*exp((v->r[i]-v->radioso)/v->aso)
		/((v->aso*v->r[i])*(1.+exp((v->r[i]-v->radioso)/v->aso))*(1.+exp((v->r[i]-v->radioso)/v->aso)));
		//misc4<<v->r[i]<<"  "<<real(potencial[i])<<"  "<<imag(potencial[i])<<endl;
	}
	funcion->wf[0]=1.e-10;
	funcion->wf[1]=(2.*(1.-0.416666667*dd*(-potencial[0]+funcion->energia))*funcion->wf[0])/
			(1.+0.083333333*dd*(-potencial[1]+funcion->energia));
	for (i=1;i<puntos-1;i++) {
		funcion->wf[i+1]=(2.*(1.-0.416666667*dd*(-potencial[i]+funcion->energia))
				*funcion->wf[i]-(1.+0.083333333*dd*(-potencial[i-1]+ funcion->energia))*funcion->wf[i-1])
				/(1.+0.083333333*dd*(-potencial[i+1]+funcion->energia));
	}
	// Calculo del desfase
	radio_1=radio_match;
	i_1=int(ceil(radio_1/delta_r))-1;
	radio_1=delta_r*(i_1 + 1.);
	fu1=funcion->wf[i_1];
	radio_2=radio_1+2.*delta_r;
	i_2=int(ceil(radio_2/delta_r))-1;
	radio_2=delta_r*(i_2+1.);
	fu2=funcion->wf[i_2];
	derivada_log=fu1/fu2;
	status1=gsl_sf_coulomb_wave_FG_e(etac,q*radio_1,funcion->l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
	//	if(status1) cout<<gsl_strerror (status1)<<"  en l="<<funcion[0].l<<",  energia="<<funcion[0].energia<<endl;
	status2=gsl_sf_coulomb_wave_FG_e(etac,q*radio_2,funcion->l,0,&F2,&Fp,&G2,&Gp,&ex1,&ex2);
	//	if(status2) cout<<gsl_strerror (status2)<<"  en l="<<funcion[0].l<<",  energia="<<funcion[0].energia<<endl;
	x=real((F1.val-F2.val*derivada_log)/(-G1.val+G2.val*derivada_log));
	y=imag((F1.val-F2.val*derivada_log)/(-G1.val+G2.val*derivada_log));
	GSL_SET_COMPLEX(&deltagsl,x,y);
	delta=(gsl_complex_arctan(deltagsl).dat[0]+I*gsl_complex_arctan(deltagsl).dat[1]); // desfase
	factor=exp(I*(delta))*(cos(delta)*F1.val+sin(delta)*G1.val)/fu1;
	//*fp<<endl<<"& Energia: "<<funcion->energia<<"  Momento angular orbital: "<<funcion->l<<"  Momento angular total: "<<funcion->j<<endl;
	for (i=0;i<puntos;i++) {
		funcion->r[i] =delta_r*(i+1.);
		funcion->wf[i]=factor*funcion->wf[i];
		if(status1 || status2) funcion->wf[i]=0.;
		gsl_sf_coulomb_wave_FG_e(etac,q*funcion->r[i],funcion[0].l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
		//*fp<<funcion->r[i]<<"  "<<real(funcion->wf[i])<<"  "<<imag(funcion->wf[i])<<endl;
	}
	delete[] potencial;
	return delta;
}

/*****************************************************************************
Funci�n de Green
 *****************************************************************************/
complejo GeneraGreenFunction(distorted_wave* funcion_regular,distorted_wave* funcion_irregular,potencial_optico *v, double q1q2,
		double masa, double radio_max,int puntos,double radio_match,double spin)
 {
	int i, i_1, i_2;
	double hbarx, dd, radio_1, radio_2, x, y, ex1, ex2, etac, q,delta_r,spinorbit;
	complejo delta,derivada_log,fu1,fu2,factor_regular;
	gsl_complex deltagsl;
	gsl_sf_result F1,G1,F2,G2,Fp,Gp;
	complejo *potencial=new complejo[puntos];
	if (funcion_regular[0].energia < 0.) Error( "Niveles intermedios energ�ticamente prohibidos!");
	delta_r=radio_max/double(puntos);
	if(radio_match>radio_max-4.*delta_r) Error("Radio de matching demasiado grande en GeneraGreenFunction");
	hbarx=HC*HC/(2.*AMU*masa);
	dd=delta_r*delta_r/hbarx;
	q=sqrt(funcion_regular[0].energia/hbarx);
	etac=q1q2*masa*E2HC*AMU/(HC*q);
	funcion_regular[0].puntos=puntos;
	funcion_regular[0].radio=radio_max;
	funcion_regular[0].j=funcion_regular[0].l+spin;
	funcion_regular[1].puntos=puntos;
	funcion_regular[1].radio=radio_max;
	funcion_regular[1].j=funcion_regular[0].l-spin;
    if(funcion_regular[1].l==0) funcion_regular[1].j=spin;

	funcion_irregular[0].puntos=puntos;
	funcion_irregular[0].radio=radio_max;
	funcion_irregular[0].j=funcion_irregular[0].l+spin;
	funcion_irregular[1].puntos=puntos;
	funcion_irregular[1].radio=radio_max;
	funcion_irregular[1].j=funcion_irregular[0].l-spin;
	if(funcion_irregular[1].l==0) funcion_irregular[1].j=spin;
	// ********************************************** C�lculo para j=l+1/2 ***********************************************************
	/* actualizacion del potencial con los t�rminos de Coulomb, centr�fugo y de spin-�rbita*/
	spinorbit =(funcion_regular[0].j*(funcion_regular[0].j+1.)-funcion_regular[0].l*(funcion_regular[0].l+1.)-spin*(spin+1.));
	for (i=0;i<puntos-1;i++) {
		if(v->r[i]>v->radio_coul) potencial[i]=v->pot[i]+E_CUADRADO*q1q2/v->r[i]+
				(funcion_regular[0].l*(funcion_regular[0].l+1.))*hbarx /(v->r[i]*v->r[i])
				-2.*spinorbit*v->Vso*exp((v->r[i]-v->radioso)/v->aso)
		/((v->aso*v->r[i])*(1.+exp((v->r[i]-v->radioso)/v->aso))*(1.+exp((v->r[i]-v->radioso)/v->aso)));
		if(v->r[i]<=v->radio_coul) potencial[i]=v->pot[i]+E_CUADRADO*q1q2*(3.-(v->r[i]/v->radio_coul)* (v->r[i]/v->radio_coul))/(2.*v->radio_coul)+
				(funcion_regular[0].l*(funcion_regular[0].l+1.))*hbarx /(v->r[i]*v->r[i])
				-2.*spinorbit*v->Vso*exp((v->r[i]-v->radioso)/v->aso)
		/((v->aso*v->r[i])*(1.+exp((v->r[i]-v->radioso)/v->aso))*(1.+exp((v->r[i]-v->radioso)/v->aso)));
//		misc1<<delta_r*(i+1.)<<"   "<<v->pot[i]<<"   "<<E_CUADRADO*q1q2/v->r[i]<<"   "<<2.*spinorbit*v->Vso*exp((v->r[i]-v->radioso)/v->aso)
//		/((v->aso*v->r[i])*(1.+exp((v->r[i]-v->radioso)/v->aso))*(1.+exp((v->r[i]-v->radioso)/v->aso)))<<endl;
//		misc1<<v->r[i]<<"   "<<real(potencial[i])<<"   "<<imag(potencial[i])<<endl;
	}
	// ********************************************* Solucion regular j=l+1/2 **********************************************************
	funcion_regular[0].wf[0]=1.e-10;
	funcion_regular[0].wf[1]=(2.*(1.-0.416666667*dd*(-potencial[0]+funcion_regular[0].energia))*funcion_regular[0].wf[0])
			/(1.+0.083333333*dd*(-potencial[1]+funcion_regular[0].energia));
	for (i=1;i<puntos-1;i++) {
		funcion_regular[0].wf[i+1]=(2.*(1.-0.416666667*dd*(-potencial[i]+funcion_regular[0].energia))
				*funcion_regular[0].wf[i]-(1.+0.083333333*dd*(-potencial[i-1]+funcion_regular[0].energia))*funcion_regular[0].wf[i-1])
				/(1.+0.083333333*dd*(-potencial[i+1]+funcion_regular[0].energia));
	}

	radio_1=radio_match;
	i_1=int(ceil(radio_1/delta_r))-1;
	radio_1=delta_r*(i_1 + 1.);
	fu1=funcion_regular[0].wf[i_1];
	radio_2=radio_1+2.*delta_r;
	i_2=int(ceil(radio_2/delta_r))-1;
	radio_2=delta_r*(i_2+1.);
	fu2=funcion_regular[0].wf[i_2];
	derivada_log=fu1/fu2;
	gsl_sf_coulomb_wave_FG_e(etac,q*radio_1,funcion_regular[0].l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
	gsl_sf_coulomb_wave_FG_e(etac,q*radio_2,funcion_regular[0].l,0,&F2,&Fp,&G2,&Gp,&ex1,&ex2);

	x=real((F1.val-F2.val*derivada_log)/(-G1.val+G2.val*derivada_log));
	y=imag((F1.val-F2.val*derivada_log)/(-G1.val+G2.val*derivada_log));
	GSL_SET_COMPLEX(&deltagsl,x,y);
	delta=(gsl_complex_arctan(deltagsl).dat[0]+I*gsl_complex_arctan(deltagsl).dat[1]); // desfase

	// ********************************************* Solucion irregular  j=l+1/2**********************************************************
	for (i =0;i<puntos;i++) {
		funcion_irregular[0].wf[i]= 0.;
	}
	gsl_sf_coulomb_wave_FG_e(etac,q*radio_max,funcion_irregular[0].l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
	gsl_sf_coulomb_wave_FG_e(etac,q*(radio_max-delta_r),funcion_irregular[0].l,0,&F2,&Fp,&G2,&Gp,&ex1,&ex2);
	funcion_irregular[0].wf[puntos-1]=I*F1.val+G1.val;
	funcion_irregular[0].wf[puntos-2]=I*F2.val+G2.val;
	for (i=puntos-2;i>0;i--) {
		funcion_irregular[0].wf[i-1]= (2.*(1.-0.416666667*dd*(-potencial[i]+funcion_irregular[0].energia))
				*funcion_irregular[0].wf[i]-(1.+0.083333333*dd*(-potencial[i+1]+funcion_irregular[0].energia))
				*funcion_irregular[0].wf[i+1])/(1.+0.083333333*dd*(-potencial[i-1]+funcion_irregular[0].energia));
	}
	gsl_sf_coulomb_wave_FG_e(etac,q*radio_1,funcion_irregular[0].l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
	factor_regular=exp(I*(delta))*(cos(delta)*F1.val+sin(delta)*G1.val)/fu1;
	for (i=0;i<puntos;i++) {
		funcion_regular[0].r[i]=delta_r*(i+1.);
		funcion_irregular[0].r[i]=delta_r*(i+1.);
		gsl_sf_coulomb_wave_FG_e(etac,q*funcion_regular[0].r[i],funcion_irregular[0].l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
		funcion_regular[0].wf[i]= factor_regular*funcion_regular[0].wf[i];
//		misc1<<funcion_regular[0].r[i]<<"  "<<real(funcion_regular[0].wf[i])
//				<<"  "<<imag(funcion_regular[0].wf[i])<<endl;
//		misc2<<funcion_regular[0].r[i]<<"  "<<real(funcion_irregular[0].wf[i])
//				<<"  "<<imag(funcion_irregular[0].wf[i])<<endl;
	}



	// ********************************************** C�lculo para j=l-1/2 ***********************************************************
	/* actualizacion del potencial con los t�rminos de Coulomb, centr�fugo y de spin-�rbita*/
	spinorbit =(funcion_regular[1].j*(funcion_regular[1].j+1.)-funcion_regular[1].l*(funcion_regular[1].l+1.)-spin*(spin+1.));
	if (funcion_regular[1].l==0) spinorbit=0.;
	for (i=0;i<puntos-1;i++) {
		if(v->r[i]>v->radio_coul) potencial[i]=v->pot[i]+E_CUADRADO*q1q2/v->r[i]+
				(funcion_regular[0].l*(funcion_regular[0].l+1.))*hbarx /(v->r[i]*v->r[i])
				-2.*spinorbit*v->Vso*exp((v->r[i]-v->radioso)/v->aso)
		/((v->aso*v->r[i])*(1.+exp((v->r[i]-v->radioso)/v->aso))*(1.+exp((v->r[i]-v->radioso)/v->aso)));
		if(v->r[i]<=v->radio_coul) potencial[i]=v->pot[i]+E_CUADRADO*q1q2*(3.-(v->r[i]/v->radio_coul)* (v->r[i]/v->radio_coul))/(2.*v->radio_coul)+
				(funcion_regular[0].l*(funcion_regular[0].l+1.))*hbarx /(v->r[i]*v->r[i])
				-2.*spinorbit*v->Vso*exp((v->r[i]-v->radioso)/v->aso)
		/((v->aso*v->r[i])*(1.+exp((v->r[i]-v->radioso)/v->aso))*(1.+exp((v->r[i]-v->radioso)/v->aso)));
//		misc2<<v->r[i]<<"  "<<real(potencial[i])<<"  "<<imag(potencial[i])<<endl;

	}
	// ********************************************* Solucion regular j=l-1/2 **********************************************************
	funcion_regular[1].wf[0]=1.e-10;
	funcion_regular[1].wf[1]=(2.*(1.-0.416666667*dd*(-potencial[0]+funcion_regular[1].energia))*funcion_regular[1].wf[0])
			/(1.+0.083333333*dd*(-potencial[1]+funcion_regular[1].energia));
	for (i=1;i<puntos-1;i++) {
		funcion_regular[1].wf[i+1]=(2.*(1.-0.416666667*dd*(-potencial[i]+funcion_regular[1].energia))
				*funcion_regular[1].wf[i]-(1.+0.083333333*dd*(-potencial[i-1]+funcion_regular[1].energia))*funcion_regular[1].wf[i-1])
				/(1.+0.083333333*dd*(-potencial[i+1]+funcion_regular[1].energia));
	}

	radio_1=radio_match;
	i_1=int(ceil(radio_1/delta_r))-1;
	radio_1=delta_r*(i_1 + 1.);
	fu1=funcion_regular[1].wf[i_1];
	radio_2=radio_1+2.*delta_r;
	i_2=int(ceil(radio_2/delta_r))-1;
	radio_2=delta_r*(i_2+1.);
	fu2=funcion_regular[1].wf[i_2];
	derivada_log=fu1/fu2;
	gsl_sf_coulomb_wave_FG_e(etac,q*radio_1,funcion_regular[1].l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
	gsl_sf_coulomb_wave_FG_e(etac,q*radio_2,funcion_regular[1].l,0,&F2,&Fp,&G2,&Gp,&ex1,&ex2);

	x=real((F1.val-F2.val*derivada_log)/(-G1.val+G2.val*derivada_log));
	y=imag((F1.val-F2.val*derivada_log)/(-G1.val+G2.val*derivada_log));
	GSL_SET_COMPLEX(&deltagsl,x,y);
	delta=(gsl_complex_arctan(deltagsl).dat[0]+I*gsl_complex_arctan(deltagsl).dat[1]); // desfase

	// ********************************************* Solucion irregular j=l-1/2 **********************************************************
	for (i =0;i<puntos;i++) {
		funcion_irregular[1].wf[i]= 0.;
	}

	gsl_sf_coulomb_wave_FG_e(etac,q*radio_max,funcion_irregular[1].l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
	gsl_sf_coulomb_wave_FG_e(etac,q*(radio_max-delta_r),funcion_irregular[1].l,0,&F2,&Fp,&G2,&Gp,&ex1,&ex2);
	funcion_irregular[1].wf[puntos-1]=I*F1.val+G1.val;
	funcion_irregular[1].wf[puntos-2]=I*F2.val+G2.val;
	for (i=puntos-2;i>0;i--) {
		funcion_irregular[1].wf[i-1]= (2.*(1.-0.416666667*dd*(-potencial[i]+funcion_irregular[1].energia))
				*funcion_irregular[1].wf[i]-(1.+0.083333333*dd*(-potencial[i+1]+funcion_irregular[1].energia))
				*funcion_irregular[1].wf[i+1])/(1.+0.083333333*dd*(-potencial[i-1]+funcion_irregular[1].energia));
	}
	gsl_sf_coulomb_wave_FG_e(etac,q*radio_1,funcion_irregular[1].l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
	factor_regular=exp(I*(delta))*(cos(delta)*F1.val+sin(delta)*G1.val)/fu1;
	for (i=0;i<puntos;i++) {
		funcion_regular[1].r[i]=delta_r*(i+1.);
		funcion_irregular[1].r[i]=delta_r*(i+1.);
		gsl_sf_coulomb_wave_FG_e(etac,q*funcion_regular[1].r[i],funcion_irregular[1].l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
		funcion_regular[1].wf[i]= factor_regular*funcion_regular[1].wf[i];
	}
	delete[] potencial;
	return delta;
}
void GeneraEstadosPI(potencial* pot,estado* st,double radio,int puntos,double cargas,parametros* parm,int ajuste,double masa
		,double* D0,double* rms)
{
	double energia,etrial,vmax,vmin;
	if(ajuste==1)
	{
		if(*(st->file)=='\0')
		{
			energia=st->energia;
			etrial=MAX_ENERGIA;
			vmax=MAX_ENERGIA;
			vmin=-MAX_ENERGIA;
			pot->V=MAX_ENERGIA;
			cout<<"energia nivel: "<<energia<<"   l: "<<st->l
					<<"   nodos: "<<st->nodos<<"   j: "<<st->j<<endl;
			while(fabs(etrial-energia)>EPSILON)
			{
//				cout<<etrial<<"  "<<energia<<endl;
				pot->V=-(vmax+vmin)/2.;
//				cout<<pot->V<<endl;
				GeneraPotencialCM(parm,pot);
//				cout<<"quillo1"<<endl;
				GeneraEstado(st,pot,radio,puntos,cargas,masa,D0,rms);
				etrial=st->energia;
//				cout<<etrial<<"  "<<pot->V<<endl;
				if(etrial>energia) vmax=-pot->V;
				if(etrial<=energia) vmin=-pot->V;
			}
		}
		else File2State(st,parm);
	}

	if(ajuste==0)
	{
		if(*(st->file)=='\0')
		{
			GeneraEstado(st,pot,radio,puntos,cargas,masa,D0,rms);
			cout<<"energia nivel: "<<st->energia<<"   l: "<<st->l
					<<"   nodos: "<<st->nodos<<"   j: "<<st->j<<endl;
		}
		else File2State(st,parm);
	}
}
/////////////////////////////////////////////////////////////////////////
//                                                                     //
//     Genera los estados de part�cula independiente de un nucleo dado //
//        y el     correspondiente potencial de campo medio            //
//                                                                     //
/////////////////////////////////////////////////////////////////////////
void GeneraEstadosContinuo(potencial_optico* pot_optico,estado* st,double radio,int puntos,double cargas,parametros* parm,double masa)
{
	double energia,hbarx,q;
	int i;
	complejo estado_inicial;
	distorted_wave* dw=new distorted_wave;
	ofstream fp1("aux.txt");
	GeneraPotencialOptico(parm,pot_optico,parm->m_A,parm->m_b);
	hbarx=HC*HC/(2.*AMU*masa);
	q=sqrt(st->energia/hbarx);
	dw->energia=st->energia;
	dw->l=st->l;
	dw->spin=parm->n_spin;
	dw->j=st->j;
	dw->radio=st->radio;
	if(*(st->file)=='\0')
	{
		GeneraDWspin(dw,pot_optico,cargas,masa,radio,puntos,parm->matching_radio,&fp1);
		st->puntos=puntos;
		cout<<"energia nivel: "<<st->energia<<"   l: "<<st->l
				<<"   nodos: "<<st->nodos<<"   j: "<<st->j<<endl;

		for (i=0;i<puntos;i++) {
			st->wf[i]=dw->wf[i]/(q*dw->r[i]);
			st->r[i]=dw->r[i];
		}
	}
	else File2State(st,parm);
}
void GeneraRemnant(potencial_optico *pot,potencial_optico *core,potencial_optico *in_pot,
		potencial_optico *in_core,double q1q2_pot,double q1q2_core,int l_pot,int l_core,double masa_pot,int masa_core)
{
	int i;
	double hbarx_pot,hbarx_core;
	hbarx_pot=HC*HC/(2.*AMU*masa_pot);
	hbarx_core=HC*HC/(2.*AMU*masa_core);
	l_pot=0;
	l_core=0;
//	misc1<<"+++++++++++++++++++++++++++++++++++++++++++"<<endl<<endl;
	for (i=0;i<in_pot->puntos;i++)
	{
		pot->r[i]=in_pot->r[i];
		core->r[i]=in_core->r[i];
		if(core->r[i]>=in_core->radio_coul) core->pot[i]=in_core->pot[i]+E_CUADRADO*q1q2_core/core->r[i]+
				(l_core*(l_core+1.))*hbarx_core /(core->r[i]*core->r[i]);
		if(core->r[i]<in_core->radio_coul) core->pot[i]=in_core->pot[i]+E_CUADRADO*q1q2_core*(3.-(core->r[i]/in_core->radio_coul)
				* (core->r[i]/in_core->radio_coul))/(2.*in_core->radio_coul)+
				(l_core*(l_core+1.))*hbarx_core /(core->r[i]*core->r[i]);
		if(pot->r[i]>=in_pot->radio_coul) pot->pot[i]=in_pot->pot[i]+E_CUADRADO*q1q2_pot/pot->r[i]+
				(l_pot*(l_pot+1.))*hbarx_pot /(pot->r[i]*pot->r[i]);
		if(pot->r[i]<in_pot->radio_coul) pot->pot[i]=in_pot->pot[i]+E_CUADRADO*q1q2_pot*(3.-(pot->r[i]/in_pot->radio_coul)
				* (pot->r[i]/in_pot->radio_coul))/(2.*in_pot->radio_coul)+
				(l_pot*(l_pot+1.))*hbarx_pot /(pot->r[i]*pot->r[i]);
//		misc1<<pot->r[i]<<"  "<<real(core->pot[i])<<"  "<<imag(core->pot[i])<<"  "<<real(pot->pot[i])<<"  "<<imag(pot->pot[i])<<endl;
	}
	core->puntos=in_core->puntos;
	pot->puntos=in_pot->puntos;
}
//////////////////////////////////////////////
//    Inicializador de tensores  complejos  //
//////////////////////////////////////////////
complejo*** tensor_cmpx(int dim1,int dim2,int dim3)
{
  int n,m,p;
  complejo*** tensor=new complejo**[dim1];
  for(n=0;n<dim1;n++)
  {
    tensor[n]=new complejo*[dim2];
    for(m=0;m<dim2;m++)
    {
      tensor[n][m]=new complejo[dim3];
    }
  }
  for(n=0;n<dim1;n++)
  {
    for(m=0;m<dim2;m++)
    {
      for(p=0;p<dim3;p++)
      {
        tensor[n][m][p]=0.;
      }
    }
  }
  return (tensor);
}
//////////////////////////////////////////////
//    Inicializador de tensores  double     //
//////////////////////////////////////////////
double*** tensor_dbl(int dim1,int dim2,int dim3)
{
  int n,m,p;
  double*** tensor=new double**[dim1];
  for(n=0;n<dim1;n++)
  {
    tensor[n]=new double*[dim2];
    for(m=0;m<dim2;m++)
    {
      tensor[n][m]=new double[dim3];
    }
  }
  for(n=0;n<dim1;n++)
  {
    for(m=0;m<dim2;m++)
    {
      for(p=0;p<dim3;p++)
      {
        tensor[n][m][p]=0.;
      }
    }
  }
  return (tensor);
}
////////////////////////////////////////////////////////
//    Inicializador de tensores de grado 4 complex     //
////////////////////////////////////////////////////////
complejo**** tensor4_cmpx(int dim1,int dim2,int dim3,int dim4)
{
  int n,m,p,l;
  complejo**** tensor4=new complejo***[dim1];
  for(n=0;n<dim1;n++)
  {
    tensor4[n]=new complejo**[dim2];
    for(m=0;m<dim2;m++)
    {
      tensor4[n][m]=new complejo*[dim3];
      for(p=0;p<dim3;p++)
      {
        tensor4[n][m][p]=new complejo[dim4];
      }
    }
  }
  for(n=0;n<dim1;n++)
  {
    for(m=0;m<dim2;m++)
    {
      for(p=0;p<dim3;p++)
      {
        for(l=0;l<dim4;l++)
        {
          tensor4[n][m][p][l]=0.;
        }
      }
    }
  }
  return (tensor4);
}
void InicializaOneTrans(struct parametros* parm)
{
	double masa_proyectil,masa_blanco;
//	parm->m_B=parm->m_A+1.;
//	parm->m_b=parm->m_a-1.;
//	if (parm->m_b<1.) Error("m_b menor que 1");
	if (!strcmp(parm->proyectil,"a")) {masa_proyectil=parm->m_a; masa_blanco=parm->m_A;}
	if (!strcmp(parm->proyectil,"A")) {masa_proyectil=parm->m_A; masa_blanco=parm->m_a;}
	if ((strcmp(parm->proyectil,"A")!=0) && ((strcmp(parm->proyectil,"a")!=0))) Error("Proyectil debe ser 'a' o 'A' ");
	parm->energia_cm=(masa_blanco/(parm->m_a+parm->m_A))*parm->energia_lab;
	if(-parm->Qvalue>parm->energia_cm) Error("Energ�a de reacci�n insuficiente");
	parm->mu_Aa=(parm->m_a*parm->m_A)/(parm->m_a+parm->m_A);
	parm->mu_Bb=(parm->m_b*parm->m_B)/(parm->m_b+parm->m_B);
	parm->k_Aa=sqrt(2.*parm->mu_Aa*AMU*parm->energia_cm)/HC;
	parm->k_Bb=sqrt(2.*parm->mu_Bb*AMU*(parm->energia_cm+parm->Qvalue))/HC;
	parm->eta=parm->Z_a*parm->Z_A*E2HC*parm->mu_Aa*AMU/(HC*parm->k_Aa);
	cout<<" ma: "<<parm->m_a<<endl;
	cout<<" mB: "<<parm->m_B<<endl;
	cout<<" mb: "<<parm->m_b<<endl;
	cout<<" masa proyectil: "<<masa_proyectil<<endl;
	cout<<" mA: "<<parm->m_A<<endl;
	cout<<" masa blanco: "<<masa_blanco<<endl;
	cout<<" energia laboratorio: "<<parm->energia_lab<<" MeV"<<endl;
	cout<<" energia CM: "<<parm->energia_cm<<" MeV"<<endl;
	cout<<" Q-value: "<<parm->Qvalue<<" MeV"<<endl;
	cout<<" masa reducida canal inicial: "<<parm->mu_Aa<<endl;
	cout<<" masa reducida canal final: "<<parm->mu_Bb<<endl;
	cout<<" momento  canal inicial: "<<parm->k_Aa<<" fm^-1"<<endl;
	cout<<" momento  canal final: "<<parm->k_Bb<<" fm^-1"<<endl;
}
//////////////////////////////////////////////////////
//  Lee estado de un archivo                       //
////////////////////////////////////////////////////
void File2State(estado *st,parametros *parm)
{
	ifstream fp;
	fp.open(st->file);
	if(!fp.is_open()) {cout<<"No se pudo abrir "<<st->file<<endl; exit(0);}
	int puntos,n;
	double pos,delta_r;
	double *r=new double[MAX_PTS];
	double *wf=new double[MAX_PTS];
	double *wf_r=new double[MAX_PTS];
	delta_r=parm->radio/double(parm->puntos);
	puntos=0;
	while(!fp.eof())
	{
		puntos++;
		fp>>r[puntos];
		fp>>wf[puntos];
		wf_r[puntos]=wf[puntos]/r[puntos];
		if(puntos>=MAX_PTS) {cout<<"N�mero de puntos en "<<st->file<<" mayor que MAX_PTS"<<endl; exit(0);}
//		misc3<<r[puntos]<<"  "<<wf[puntos]<<"  "<<wf_r[puntos]<<endl;
	}

	for(n=0;n<parm->puntos;n++)
	{
		pos=delta_r*(n+1);
		st->r[n]=pos;
		if(pos<=r[puntos-1]) st->wf[n]=interpola_dbl(wf_r,r,pos,puntos-1);
		else st->wf[n]=0.;
//		misc2<<pos<<"  "<<st->wf[n]*st->r[n]<<"  "<<st->wf[n]<<endl;
	}
	st->puntos=parm->puntos;
	st->radio=parm->radio;
	delete[] r;
	delete[] wf;
	delete[] wf_r;
}
void GeneraEstado(estado *st,potencial *potencial, double radio_max,int puntos,double q1q2,double masa, double* D0
		,double* rms) {
	int ND,i;
	double *vv=new double[puntos],*sx=new double[puntos],*vs=new double[puntos],*v=new double[puntos];
	double hbarx,dd,Wlim,centr,Etrial,Emax,Emin,ls,delta_r,radio;
	delta_r=radio_max/double(puntos);
	Emin=MIN_ENERGIA;
	Emax=MAX_ENERGIA;
	hbarx=HC*HC/(2.*AMU*masa);
	dd=delta_r*delta_r/hbarx;
	Wlim=1.e-13;
	centr=(st->l*(st->l+1.))*hbarx;
	ls=(st->j*(st->j+1.)-st->l*(st->l+1.)-0.75);
	st->puntos=puntos;
	// a�ade los t�rminos Coulomb y spin-�rbita *********************************************************************
	if(!strcmp(potencial->tipo,"ws"))
	{
		for (i=0;i<puntos;i++) {
			st->r[i]=delta_r*(i+1);
			if(st->r[i]>potencial->RV) v[i]=potencial->pot[i]+E_CUADRADO*q1q2/st->r[i];
			if(potencial->r[i]<=potencial->RV) v[i]=potencial->pot[i]
								+E_CUADRADO*q1q2*(3.-(st->r[i]/potencial->RV)* (st->r[i]/potencial->RV))/(2.*potencial->RV);
			vs[i]=-2.*(potencial->VSO)*exp((st->r[i]-potencial->RSO)/potencial->aSO)/
					(st->r[i]*potencial->aSO*(1.+exp((st->r[i]-potencial->RSO)
							/potencial->aSO))*(1.+exp((st->r[i]-potencial->RSO)/potencial->aSO)));
			potencial->pot[i]=v[i];
//			misc3<<st->r[i]<<"  "<<v[i]<<"  "<<ls*vs[i]<<"  "<<(centr)/(st->r[i]*st->r[i])<<endl;
		}
	}
	if(!strcmp(potencial->tipo,"tang"))
	{
		if(delta_r>potencial->rhc/2.) Error("Paso de intagraci�n demasiado grande para el potencial Tang-Herndon");
		for (i=0;i<puntos;i++) {
			st->r[i]=delta_r*(i+1);
			if(st->r[i]>potencial->rhc) v[i]=-potencial->V*exp(-potencial->k*(st->r[i]-potencial->rhc));
			if (st->r[i]<=potencial->rhc) v[i]=0.;
		}
	}
	while (fabs(-Emin+Emax)>Wlim) {
		Etrial=(Emin+ Emax)/2.0;
		st->wf[0]=pow(delta_r,(st->l+1));
		st->wf[1]=pow(2.0*delta_r,(st->l+1));
		if(!strcmp(potencial->tipo,"tang")) {st->wf[0]=0.;st->wf[1]=0.;}
		ND=0;
		for (i=1;i<puntos-1;i++) {
			vv[i]=v[i]+ls*vs[i]+(centr)/(st->r[i]*st->r[i]);
			sx[i]=dd*(vv[i]-Etrial);
			st->wf[i+1]=(2.0+sx[i])*st->wf[i]-st->wf[i-1];
			if(!strcmp(potencial->tipo,"tang")  && st->r[i]<=potencial->rhc) {st->wf[i]=0.;st->wf[i+1]=1.e-4;}
			if (real(st->wf[i+1]*st->wf[i])<0.) ND=ND+1;
		}
		if (ND>st->nodos) Emax=Etrial;
		if (ND<=st->nodos) Emin=Etrial;
	}
	st->energia=Etrial;
	radio=0.;
	for (i=0;i<puntos;i++) {
		st->wf[i]=st->wf[i]/st->r[i];
		if (i>0)
			if((abs(st->wf[i])>abs(st->wf[i-1]))&&(abs(st->r[i])>4.*(potencial->RV))&&(st->energia<0.))
				st->wf[i]=st->wf[i-1];
	}
	st->radio=radio_max;
	//  Normalizacion
	Normaliza(st,st,radio_max,puntos,'s');
	*D0=3.54490770181103*VertexD0(st,potencial,radio_max,puntos,rms);
//	cout<<*D0<<"  "<<*rms<<endl;
	delete[] vv;
	delete[] v;
	delete[] sx;
	delete[] vs;
}
//////////////////////////////////////////////////////
//                                                  //
//             Desfase Coulombiano                  //
//////////////////////////////////////////////////////

double deltac(int l,double etac)
{
  int neuler=5000;
  int n,i;
  double euler=0.5772156649;
  double constante=0.;
  double tr=PI/180.;
  double dif,dumb;
  double* delc=new double [neuler+1];
  for(n=1; n<l+1; n++)
  {
    constante+=1./double(n);
  }
  delc[0]=etac/(1.+double(l))-atan(etac/(1.+double(l)));
  for(n=1; n<neuler; n++)
  {
    delc[n]=delc[n-1]+etac/(1.+double(l)+double(n))
            -atan(etac/(1.+double(l)+double(n)));
    dif=delc[n]-delc[n-1];
    if (abs(dif)<1.e-4) break;
  }
  dumb=delc[n]+etac*(constante-euler);
  delete[] delc;
  return dumb;
}
/*****************************************************************************
Escribe el potencial optico
 *****************************************************************************/
void EscribePotencialOptico(int puntos,potencial_optico* pot,int numero_potenciales,struct parametros *parm)
{
	ofstream fp("potenciales_opticos.txt");
	int n,i;
	for(i=0;i<puntos;i++){
		fp<<pot[0].r[i]<<"  ";
		for(n=0;n<numero_potenciales;n++){
			fp<<real(pot[n].pot[i])<<"  "<<imag(pot[n].pot[i])<<"  ";
		}
		fp<<endl;
	}
	fp.close();
}
/*****************************************************************************
Calcula el parametro D0 de zero range
 *****************************************************************************/
double VertexD0(estado* st,potencial* pot, double radio, int pts, double* rms) {
	int regla_r, nr, indice;
	double ar, br, norma, r,radio_medio,potint;
	parametros_integral *dim=new parametros_integral;
	double step=radio/double(pts);
	regla_r=60;
	complejo sum=0.;
	complejo wfint;
	dim->a=0;
	dim->b=radio;
	dim->num_puntos=regla_r;
	GaussLegendre(dim->puntos,dim->pesos,dim->num_puntos);
	radio_medio=0.;
	for (nr=0;nr<dim->num_puntos;nr++) {
		r=dim->a+(dim->b-dim->a)*(dim->puntos[nr]+1.)/2.;
		wfint=interpola_cmpx(st->wf,st->r,r,st->puntos);
		potint=interpola_dbl(pot->pot,pot->r,r,pot->puntos);
		sum+=wfint*potint*r*r*dim->pesos[nr];
		radio_medio+=abs(wfint*wfint)*r*r*r*r*dim->pesos[nr];
	}
	norma=abs(sum)*(dim->b-dim->a)/2.;
	*rms=sqrt(radio_medio*(dim->b-dim->a)/2.);
	delete[] dim;
	return norma;
}
void AddCoulomb(potencial_optico* v,double q1q2)
{
	int i;
	for (i=0;i<v->puntos-1;i++) {
			if(v->r[i]>=v->radio_coul) v->pot[i]+=E_CUADRADO*q1q2/v->r[i];
			if(v->r[i]<v->radio_coul) v->pot[i]+=E_CUADRADO*q1q2*(3.-(v->r[i]/v->radio_coul)* (v->r[i]/v->radio_coul))/(2.*v->radio_coul);
			//misc4<<v->r[i]<<"  "<<real(potencial[i])<<"  "<<imag(potencial[i])<<endl;
		}
}
