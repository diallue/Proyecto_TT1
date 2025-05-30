#include "..\include\global.hpp"

Matrix eopdata;
Matrix Cnm;
Matrix Snm;
Param AuxParam;
Matrix PC;
Matrix obs;

void eop19620101(int c) {
	eopdata = zeros(13, c);
	
	FILE *fid = fopen("../data/eop19620101.txt","r");
	
	if (fid == NULL) {
		printf("Fail open eop19620101.txt file\n");
		exit(EXIT_FAILURE);
	}
	
	for (int j = 1; j <= c; j++) {
		fscanf(fid,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &eopdata(1, j), 
		&eopdata(2, j),&eopdata(3, j),&eopdata(4, j),&eopdata(5, j),&eopdata(6, j),&eopdata(7, j),
		&eopdata(8, j),&eopdata(9, j),&eopdata(10, j),&eopdata(11, j),&eopdata(12, j),&eopdata(13, j));
	}
	
	fclose(fid);
}

void GGM03S(int n) {
	Cnm=zeros(n,n);
	Snm=zeros(n,n);
	
	FILE *fid = fopen("../data/GGM03S.txt","r");
	
	if(fid==NULL){
		cout << "Fail open GGM03S.txt file \n";
		exit(EXIT_FAILURE);
	}
	
	double aux;
	for(int i=1;i<=n;i++){
		for (int j=1;j<=i;j++){
			fscanf(fid,"%lf %lf %lf %lf %lf %lf",
				   &aux,&aux,&Cnm(i,j),&Snm(i,j),&aux,&aux);
		}
	}
	
	fclose(fid);
}

void DE430Coeff(int row,int column) {
	PC=zeros(row,column);
	
	FILE *fid = fopen("../data/DE430Coeff.txt","r");
	
	if(fid==NULL){
		cout << "Fail open DE430Coeff.txt file \n";
		exit(EXIT_FAILURE);
	}
	
	for(int i=1;i<=row;i++){
		for (int j=1;j<=column;j++){
			fscanf(fid,"%lf",
				&PC(i,j));
		}
	}
	
	fclose(fid);
}

void AuxParamLoad() {
	AuxParam.Mjd_UTC = 4.974611635416653e+04;
	AuxParam.Mjd_TT = 4.974611706231468e+04;
	AuxParam.n = 20;
	AuxParam.m = 20;
	AuxParam.sun = 1;
	AuxParam.moon = 1;
	AuxParam.planets = 1;
}

void GEOS3(int nobs) {
    obs = Matrix(nobs, 4);
    std::ifstream geos3_file("../data/GEOS3.txt");
    if (!geos3_file) {
        std::cerr << "Error: No se pudo abrir GEOS3.txt\n";
        exit(EXIT_FAILURE);
    }
    std::string line;
    for (int i = 1; i <= nobs; ++i) {
        if (!std::getline(geos3_file, line) || line.empty()) {
            std::cerr << "Error: No se pudieron leer suficientes datos de GEOS3.txt en la línea " << i << "\n";
            geos3_file.close();
            exit(EXIT_FAILURE);
        }
        try {
            int Y = std::stoi(line.substr(0, 4));
            int M = std::stoi(line.substr(5, 2));
            int D = std::stoi(line.substr(8, 2));
            int h = std::stoi(line.substr(12, 2));
            int m = std::stoi(line.substr(15, 2));
            double s = std::stod(line.substr(18, 6));
            double az = std::stod(line.substr(25, 8));
            double el = std::stod(line.substr(35, 7));
            double Dist = std::stod(line.substr(44, 10));
            obs(i, 1) = Mjday(Y, M, D, h, m, s);
            obs(i, 2) = RAD * az;
            obs(i, 3) = RAD * el;
            obs(i, 4) = 1e3 * Dist;
        } catch (const std::exception& e) {
            std::cerr << "Error: Formato inválido en GEOS3.txt en la línea " << i << ": " << e.what() << "\n";
            geos3_file.close();
            exit(EXIT_FAILURE);
        }
    }
    geos3_file.close();
}