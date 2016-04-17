#ifndef PROCESS_LIB_RICHARDSMATERIALPROPERTY_H_
#define PROCESS_LIB_RICHARDSMATERIALPROPERTY_H_

namespace ProcessLib
{

namespace RichardsFlow {


	static double getKrelbySw_van(const double Sw/*wetting saturation*/, size_t idx_phase)
    {
        double k_rw = 0.0;
		double n = 1.8;
		double m = 1.0 - 1.0 / n;
		double S_gr = 0.0;
		double S_lr = 0.0;
		double S_le = (Sw - S_lr) / (1 - S_lr - S_gr);
		k_rw = pow(S_le, 0.5)*pow((1 - pow((1 - pow(S_le, 1 / m)), m)), 2);
        return k_rw;
    }

	/**
    * return the water saturation value by capillary pressure
	* standard van - Genucheten 
	* without any regularization
	* becareful any non-physical pc
    */
	static double getSwbyPc_van(double Pc)
    {	
		/// here all the vG parameter should delivered from input file ....
        double Sw = 0.0;
		double n = 1.8;
        double m = 1.0 - 1.0 / n;
		double P0 = 0.84e+4;
		double S_gr = 0.0;
		double S_lr = 0.0;
		Sw = 1-(1 - (((1 - S_gr - S_lr) / pow((pow((Pc / P0), n) + 1), m)) + S_lr));
        return Sw; 
    }

	static double getdSwdPc_van (double Pc)
    {
        double dSwdPc = 0.0;
		double n = 1.8;
		double m = 1.0 - 1.0 / n;
		double P0 = 0.84e+4;
		double S_gr = 0.0;
		double S_lr = 0.0;
		dSwdPc = -m*n*
			(1 - S_gr - S_lr)*(1 / P0)*(pow((Pc / P0), (n - 1)))*pow(((pow((Pc / P0), n)) + 1), (-m - 1));      
        return dSwdPc; 
    }

}  // Richards_materialproperty

}  // ProcessLib


#endif  // PROCESS_LIB_RICHARDSMATERIALPROPERTY_H_