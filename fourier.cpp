#define _USE_MATH_DEFINES
#include <iostream>
#include "Wave.cpp"
#include "Wave.hpp"
#include <math.h>
#include <fstream>
#include <vector>
#include <cfloat>
#include <time.h>

void add_note(double f_note, double fe, int nb_e, double *data_real, int division=1){

	//double interval = 1/fe; // interval entre les echantillons
	long int k;
	double alpha = 2. * M_PI * f_note/fe; // constante

	for(k = 0; k < nb_e; k++){
		data_real[k] += sin((double) k * alpha)/division;
		//alpha*=0.999;
	}
}

// make note gamme chromatique
void makeNoteGc(double f_note, double fe, int nb_e, double *data_real, int start_index){

	//double interval = 1/fe; // interval entre les echantillons
	long int k;
	double alpha = 2. * M_PI * f_note/fe; // constante

	for(k = start_index; k < start_index+nb_e; k++){
		data_real[k] = sin((double) k * alpha);
		//alpha*=0.999;
	}
}

double char2double(unsigned char val){
	double real;
	real = (double) val/127.5;
	return real-1.0;
}

unsigned char double2char(double value){
	/*
	val = (val + 1.0) * 127.5;
	if(val < 0.0) val = 0.0;
	if(val > 255.0) val = 255.0;
	return (unsigned char) floor (val);*/

	if(value < -1.0){return (unsigned char)0;}
    if(value > 1.0) {return (unsigned char)255;}
    return (unsigned char)floor((value+1.0)*127.5);
}

// transformee de fourier discrete
void TFD(double *signal, double *real_part, double *imaginary_part, int nb_e){
	// signal reste le meme
	// real_part et imaginary_part representent le nouveau signal séparer en deux parties car nre complex

	double cons = 2.0 * M_PI * (1.0/nb_e); // pour gagner du temps de calcul

	int k, n;
	for(k = 0; k < nb_e; k++){
		double cons_k = cons*k; // pour gagner du temps de calcul
		imaginary_part[k] = 0;
		real_part[k] = 0;

		for(n = 0; n < nb_e; n++){
			imaginary_part[k] += signal[n] * sin(cons_k*n);
			real_part[k] += signal[n] * cos(cons_k*n);
		}

		if(k%1000 == 0){
			std::cout << "progression: " << (k/1000) << "/" << (nb_e/1000) << std::endl;
		}

		imaginary_part[k] *= (-1.0);
	}
}

// transformee de fourier inverse
void TFI(double *real_part, double *imaginary_part, double *a, double *b, int nb_e){
	
	double cons1 = 1.0/(double)nb_e;
	double cons2 = 2.0 * M_PI * cons1;

	for(long int n = 0; n < nb_e; n++){
		//imaginary_part[n] = 0.0;
		real_part[n] = 0.0;

		double cons3 = cons2*n;

		for(long int k = 0; k < nb_e; k++){
			real_part[n] += a[k]*cos(cons3*k) - b[k]*sin(cons3*k);
			//imaginary_part[n] += a[k]*sin(cons3*k) + b[k]*cos(cons3*k);
		}
		if(n%1000 == 0){
			std::cout << "progression: " << (n/1000) << "/" << (nb_e/1000) << std::endl;
		}

		//imaginary_part[n] *= cons1;
		real_part[n] *= cons1;
	}
}

void filter_passe_bas(double *real_part, double *imaginary_part, double *real_part_filter, double *imaginary_part_filter, double fc, int nb_e, double to, double fe){

	int fc_index = (int) (fc/fe*nb_e);

	for(int i = 0; i < fc_index; i++){
		real_part_filter[i] = real_part[i];
        imaginary_part_filter[i] = imaginary_part[i];
	}
	for(int j = fc_index; j < nb_e; j++){
		real_part_filter[j] = 0.0;
        imaginary_part_filter[j] = 0.0;
	}

	for(int k = nb_e - fc_index; k < nb_e; k++){
        real_part_filter[k] = real_part[k];
        imaginary_part_filter[k] = imaginary_part[k];
    }

	return ;
}

void filter_passe_haut(double *real_part, double *imaginary_part, double *real_part_filter, double *imaginary_part_filter, double fc, int nb_e, double to, double fe){

	int fc_index = (int) (fc/fe*nb_e);

	for(int i = 0; i < fc_index; i++){
        real_part_filter[i] = 0.0;
        imaginary_part_filter[i] = 0.0;
	}
	for(int j = fc_index; j < nb_e; j++){
		real_part_filter[j] = real_part[j];
        imaginary_part_filter[j] = imaginary_part[j];
	}

	for(int k = nb_e - fc_index; k < nb_e; k++){
        real_part_filter[k] = 0.0;
        imaginary_part_filter[k] = 0.0;
    }

	return;
}
void filter_coupe_bande(double *real_part, double *imaginary_part, double *real_part_filter, double *imaginary_part_filter, double fc, double fc2, int nb_e, double to, double fe){

	int fc_index = (int) (fc/fe*nb_e);
	int fc2_index = (int) (fc2/fe*nb_e);

	for(int i = 0; i < nb_e; i++){
		if (((i>=fc_index)&&(i<=fc2_index))||((i>=nb_e-fc2_index)&&(i<=nb_e-fc_index))) {
			real_part_filter[i] = 0.0;
			imaginary_part_filter[i] = 0.0;
		} else {
			
			real_part_filter[i] = real_part[i];
			imaginary_part_filter[i] = imaginary_part[i];
		}
	}
	

	return ;
}

// filtre de butterworth d'ordre 3
void filter_butterworth(double *signal, double *signal_filt, int nb_e, double fe, double fc){

	double alpha = M_PI*fc/fe;

	double A, B, C , D;
    A = (1 + 2*alpha + 2*pow(alpha,2) + pow(alpha,3));
    B = (-3 - 2*alpha + 2*pow(alpha,2) + 3*pow(alpha,3));
    C = (3 - 2*alpha - 2*pow(alpha,2) + 3*pow(alpha,3));
    D = (-1 + 2*alpha - 2*pow(alpha,2) + pow(alpha,3));

    double a[4], b[4];
    b[0] = pow(alpha,3);
    b[3] = b[0];
    b[1] = 3.0 * b[0];
    b[2] = b[1];
    a[0] = 0.0;
    a[1] = (- B/A);
    a[2] = (- C/A);
    a[3] = (- D/A);

	int n, k ,m;
	for(n = 0; n < nb_e; n++){
		signal_filt[n] = 0.0;

		for(k = 0; k < 4; k++){
			m = n-k;

			if(m>=0){
				signal_filt[n] += b[k] * signal[m];
                signal_filt[n] += a[k] * signal_filt[m];
			}
		}
	}
}

void export_wav_file(double data[], // Tableau de don�es lorsque l'on est sur des donn�es 8 bits
           long int data_nb,     // Nombre de donn�es
           short channels_nb,    // Nombre de canaux (1 pour mono ou 2 pour st�r�o)
           int sampling_freq,   // frequence d'echantillonage
		   string name   		// nom de l'onde analysee pour le fichier
		   ) {	

	unsigned char *data8;
	data8 = new unsigned char[data_nb];
	for(long int i = 0; i < data_nb; i++){
		data8[i] = double2char(data[i]);
	}
	Wave domisol = Wave(data8, data_nb, channels_nb, sampling_freq);
	//string domisol_filename = string(name) + freq_echant + string("_")+ duration + string(".wav");
	//const char* domisol_file = domisol_filename.c_str();
	domisol.write(name.c_str());
	std::cout << "fichier " << name << " created." << std::endl;

}

double *normalized(long int nb_e,double real_part[], double imaginary_part[]){
	double *data_norms;
	data_norms = new double[nb_e];
	double *fourier_data;
	fourier_data= new double[nb_e];

	double max_norm = -DBL_MAX;
	for(long int i = 0; i < nb_e; i++){
		data_norms[i] = sqrt(real_part[i]*real_part[i] + imaginary_part[i]*imaginary_part[i]);
		if(data_norms[i] > max_norm){
			max_norm = data_norms[i];
		}
	}
	for(long int i = 0; i < nb_e; i++){
		fourier_data[i] = data_norms[i]/max_norm;
	}
	delete [] data_norms;
	return fourier_data;
}
double *make_chroma(double to, double fe, double gamme[],int gamme_size){
	//std::vector<double> notes_freq = {440., 494., 523., 587., 659., 698., 784.};
	//double fe = 11000.;
	//double to = 2.;
	//double fe = 44100.;
	//double to = 7.;
	//int gamme_size = gamme.size()/sizeof(gamme[0]);
	std::cout << "size of gamme = " << gamme_size << " note(s)" << std::endl;
	long int nb_e = (long int) floor(to*fe); 
	long int nb_e_n = (long int) floor(nb_e/gamme_size); 
	
	double *data_chroma;
	data_chroma = new double[nb_e];

	for(int i=0; i < gamme_size; i++){
		//double freq_note = gamme[i];
		std::cout << "Note #" << i << std::endl;
		std::cout << "Appel de la méthode MakeNoteGc" << std::endl;
		makeNoteGc(gamme[i], fe, nb_e_n, data_chroma, i*nb_e_n);
	}
	return data_chroma;
}

double *make_accord(double to, double fe, double accord[], int accord_size){
	short channels_nb = 1;
	
	std::cout << "\nsize of accord = " << accord_size << " note(s)" << std::endl;
	
	const long int nb_e = (long int) floor(to*fe);
	
	double *data_accord;
	data_accord = new double[nb_e];

	for (int z=0;z<nb_e;z++){
		data_accord[z] = 0.;
	}
	for (int n=0;n < accord_size;n++){
		add_note(accord[n], fe, nb_e, data_accord, accord_size);
	}

	return data_accord;
}


void make_accordWithTFD(double to, double fe, double accord[], int accord_size){
	short channels_nb = 1;
	
	std::cout << "\nsize of accord = " << accord_size << " note(s)" << std::endl;
	
	const long int nb_e = (long int) floor(to*fe);

	double *TFD_accord_real;
	TFD_accord_real = new double[nb_e];

	double *TFD_accord_imaginary;
	TFD_accord_imaginary = new double[nb_e];

	double *one_note;
	one_note = new double [1];

	double *data_one_note;
	data_one_note = new double[nb_e];

	double *TFD_one_note_real;
	TFD_one_note_real = new double[nb_e];
	double *TFD_one_note_imaginary;
	TFD_one_note_imaginary = new double[nb_e];

	double *real_sum;
	real_sum = new double[nb_e];
	double *imaginary_sum;
	imaginary_sum = new double[nb_e];

	for (int i=0;i < accord_size;i++){
		one_note[0] = accord[i];
		data_one_note = make_chroma(to, fe , one_note, 1);
		TFD(data_one_note, TFD_one_note_real, TFD_one_note_imaginary, nb_e);

		for(int j = 0; j < nb_e; j++){
			real_sum[j] += TFD_one_note_real[j];
			imaginary_sum[j] += TFD_one_note_imaginary[j];
		}
	}

	for(int k=0; k< nb_e; k++){
		TFD_accord_real[k] = real_sum[k]/accord_size;
		TFD_accord_imaginary[k] = imaginary_sum[k]/accord_size;
	}

	// Proceed with TFI
	double *real_part_i;
	real_part_i = new double[nb_e];
	double *imaginary_part_i;
	imaginary_part_i = new double[nb_e];

	TFI(real_part_i, imaginary_part_i, TFD_accord_real, TFD_accord_imaginary, nb_e);
	
	// Export to Wav file
	string filename = string("Accord_from_TFD") + string(".wav");
	export_wav_file(real_part_i, nb_e, channels_nb, fe, filename);
}

double *read_piano(long int &nb_e, double &to, double &fe, string filename){
	// ----------------------------------------------------------------
	// Créer un tableau avec le fichier gamme piano
	// ----------------------------------------------------------------
	Wave vraie_gamme = Wave();
	//string filename = "GammePiano.wav";
	//char* vraie_gamme_name = filename.c_str();
	vraie_gamme.read(filename.c_str());
	unsigned char* data8;
	//int nb_e;
	vraie_gamme.getData8(&data8, &nb_e);

	double *data_read;
	data_read = new double[nb_e];
	for(long int i = 0; i < nb_e; i++){
		data_read[i] = char2double(data8[i]);
	}
	//double fe =  44100.;

	int fe_int = vraie_gamme.sampling_freq; 
    double fe_inverse = 1.0 / (double)(fe_int);
     to = fe_inverse * nb_e;
     fe = (double) fe_int;

	return data_read;
}

int main(int argc, char** argv){
	//
	// Arguments
	// 		- Type of Wave A: Accord, C: Chroma, P: Piano 
	//		- Duration in sec
	//		- Sampling Frequency
	//		- Filter 0: None, 1: Low-Pass, 2: High-Pass, 3: Butter, 4: Pass-Band, 5: Cut-Band
	//		- Cutoff Frequency of filter
	//		- Second cutoff Frequency of filter (for pass-band and cut-band)
	//

	// default values
	double to = 1.;
	double fe = 22050.;
	double fc = 350.;
	double fc2 = 700.;

	int start_time = time(NULL);

	string filename = "file";
	string duration = "1sec";
	string frequence_echantillonage = "22050Hz";
	string filtered = "None";
	string cutoff_freq = "--";
	string cutoff_freq2 = "";
	string wave_type = "A";
	string type_name = "Accord_";

	std::cout << "Fourier v1.0 \n" << std::endl;
	
	// Arguments
	cout << "You have entered " << argc
         << " arguments:" << "\n";
  
    for (int i = 0; i < argc; ++i)
        cout << argv[i] << "\n";
  
    
	if (argc!=7){
		std::cout << "Invalid number of arguments, using default" << std::endl; 
		std::cout << "Usage: fourier [type] [duration] [sampling_frequency] [1 for Low, 2 for High, 3 for Butter, 4 for Band, 5 for Cutband, 0 for None] [cutoff frequency] [cutoff frequency 2]" << std::endl; 
	}
	else {
		wave_type = argv[1];
		to = stod(argv[2]);
		fe = stod(argv[3]);
		int filter = atoi(argv[4]);
		std::cout << "Filter = " << filter << std::endl; 
		switch(filter) {
		case 3:
			// code block
			filtered = "Butter";
			break;

		case 2:
			// code block
			filtered = "High";
			break;
		
		case 1:
			// code block
			filtered = "Low";
			break;
		case 4:
			// code block
			filtered = "Band";
			break;
		case 5:
			// code block
			filtered = "Cutband";
			break;
		default:
			// code block
			filtered = "None";
		}	
		fc = stod(argv[5]);	// frequence coupure 1
		fc2 = stod(argv[6]);	// frequence coupure 2 for pass-band

		duration = string(argv[2])+string("sec");
		frequence_echantillonage = string(argv[3])+string("Hz");
		cutoff_freq = string(argv[(5)])+string("Hz");
		cutoff_freq2 = string(argv[(6)])+string("Hz");
		
		
	}
	std::cout << "Duration: " + duration + " - Echantillonage: " + frequence_echantillonage + " - Filter: " + filtered + " - Cutoff: " + cutoff_freq << std::endl; 


	short channels_nb = 1;
	
	long int nb_e = (long int) floor(to*fe);
	double *data_total;
	data_total = new double[nb_e];

	// set of notes of the accord (do mi sol la)
	//double notes_freq_accord[] = {261.,329.,392.,440.};
	// do mi sol
	//double notes_freq_accord[] = {261.,329., 392.};
	// la do mi
	double notes_freq_accord[] = {440., 523., 659.};
	// single note accord LA
	//double notes_freq_accord[] = {523.};  // LA = 440, DO= 523
	// Gamme Chromatique en La
	double notes_freq_chroma[] = {440., 494., 523., 587., 659., 698., 784.};
	
	
	
	if (wave_type=="A"){
		// make the accord with the requested notes
		type_name = "Accord_";
		int n_notes = sizeof(notes_freq_accord)/sizeof(notes_freq_accord[0]);
		data_total = make_accord(to, fe, notes_freq_accord, n_notes);
	} else if (wave_type=="C"){
		// makes the gamme chroma
		type_name = "Chroma_";
		int n_notes = sizeof(notes_freq_chroma)/sizeof(notes_freq_chroma[0]);
		data_total = make_chroma(to, fe, notes_freq_chroma, n_notes);
	} else if (wave_type=="P"){
		data_total = read_piano(nb_e,to,fe,string("GammePiano.wav"));
		type_name = "Piano_";
		duration = to_string(int(to))+string("sec");
		frequence_echantillonage = to_string(int(fe))+string("Hz");
	} else if (wave_type=="C2"){
		int n_notes = sizeof(notes_freq_accord)/sizeof(notes_freq_accord[0]);
		make_accordWithTFD(to, fe, notes_freq_accord, n_notes);
		int end_time = time(NULL);
		int elapsed_time = end_time - start_time;
		std::cout << "Temps de calcul: "  << elapsed_time << " sec" << std::endl;
		return 0;
	}
	// export to wav file
	filename = type_name + string("Filter_") + filtered + string("_cutoff_")+ cutoff_freq + string("_")+ cutoff_freq2 + string("_echant_")+ frequence_echantillonage + string("_")+ duration + string(".wav");
	export_wav_file(data_total, nb_e, channels_nb, fe, filename);


	// Proceed with TFD on data_total
	double *real_part;
	real_part = new double[nb_e];
	double *imaginary_part;
	imaginary_part = new double[nb_e];
	std::cout << "Appel de la méthode TFD" << std::endl;
	// Fourier Transform into real_part and imaginary_part
	TFD(data_total, real_part, imaginary_part, nb_e);

	
	// Export to wav file
	filename = string("TFD_") + type_name + filtered + string("_cutoff_")+ cutoff_freq + string("_")+ cutoff_freq2 + string("_echant_")+ frequence_echantillonage + string("_")+ duration + string(".wav");
	export_wav_file(normalized(nb_e, real_part, imaginary_part), nb_e, channels_nb, fe,filename);

	// filtre non destructif -> data d'origine real_part et imaginary_part conservees
	double *real_part_filter;
	real_part_filter = new double[nb_e];
	double *imaginary_part_filter;
	imaginary_part_filter = new double[nb_e];
	double *real_part_filter2;
	real_part_filter2 = new double[nb_e];
	double *imaginary_part_filter2;
	imaginary_part_filter2 = new double[nb_e];
	double * data_total_filtered;
	data_total_filtered = new double[nb_e];
	
	std::cout << "Filtering..." << std::endl;
	if (filtered=="Low"){
		filter_passe_bas(real_part, imaginary_part, real_part_filter, imaginary_part_filter, fc, nb_e, to, fe);
	} else if (filtered=="High") {
		filter_passe_haut(real_part, imaginary_part, real_part_filter, imaginary_part_filter, fc, nb_e, to, fe);
	} else if (filtered=="Band") {
		// first pass bas with highest cutoff freq, then pass haut with lowest  cutoff freq
		filter_passe_bas(real_part, imaginary_part, real_part_filter2, imaginary_part_filter2, fc2, nb_e, to, fe);	
		filter_passe_haut(real_part_filter2, imaginary_part_filter2, real_part_filter, imaginary_part_filter, fc, nb_e, to, fe);
	}
	else if (filtered=="Cutband") {		
		filter_coupe_bande(real_part, imaginary_part, real_part_filter, imaginary_part_filter, fc, fc2, nb_e, to, fe);
	}
	else if (filtered == "Butter") {
		filter_butterworth(data_total, data_total_filtered, nb_e, fe, fc);
	}else {
		// passe-bas with cutoff at fe = passe-TOUT fc = fe
		filter_passe_bas(real_part, imaginary_part, real_part_filter, imaginary_part_filter, fe, nb_e, to, fe);
	}
	// 
	if(filtered == "Butter"){
		filename = type_name + filtered + string("_cutoff_")+ cutoff_freq + string("_echant_")+ frequence_echantillonage + string("_")+ duration + string(".wav");
		export_wav_file(data_total_filtered, nb_e, channels_nb, fe,filename);
		std::cout << "Appel de la méthode TFD" << std::endl;
		// Fourier Transform into real_part and imaginary_part of the Butter filtered data
		TFD(data_total_filtered, real_part, imaginary_part, nb_e);
		// Export to wav file
		filename = string("Filtered_TFD_") + type_name + filtered + string("_cutoff_")+ cutoff_freq + string("_echant_")+ frequence_echantillonage + string("_")+ duration + string(".wav");
		export_wav_file(normalized(nb_e, real_part, imaginary_part), nb_e, channels_nb, fe,filename);

	}else {
		// Export to wav file
		filename = string("Filtered_TFD_") + type_name + filtered + string("_cutoff_")+ cutoff_freq + string("_")+ cutoff_freq2  + string("_echant_")+ frequence_echantillonage + string("_")+ duration + string(".wav");
		export_wav_file(normalized(nb_e, real_part_filter, imaginary_part_filter), nb_e, channels_nb, fe,filename);

		// Proceed with TFI
		double *real_part_i;
		real_part_i = new double[nb_e];
		double *imaginary_part_i;
		imaginary_part_i = new double[nb_e];
		std::cout << "Appel de la méthode TFI" << std::endl;
		TFI(real_part_i, imaginary_part_i, real_part_filter, imaginary_part_filter, nb_e);
		
		// Export to Wav file
		filename = string("TFI_") + type_name + filtered + string("_cutoff_")+ cutoff_freq + string("_")+ cutoff_freq2 + string("_echant_")+ frequence_echantillonage + string("_")+ duration + string(".wav");
		export_wav_file(real_part_i, nb_e, channels_nb, fe,filename);

		delete[] real_part_i;
		delete[] imaginary_part_i;
	}


	delete[] data_total;
	delete[] real_part_filter;
	delete[] imaginary_part_filter;
	delete[] real_part_filter2;
	delete[] imaginary_part_filter2;
	delete[] data_total_filtered;
	//delete[] fourier_data8_i;
	//delete[] data_norms;

	int end_time = time(NULL);
	int elapsed_time = end_time - start_time;
	std::cout << "Temps de calcul: "  << elapsed_time << " sec" << std::endl;

	return 0;
}