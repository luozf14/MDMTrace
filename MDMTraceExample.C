#include "MDMTrace.h"
#include <iostream>
#include <fstream>
#include "json.h"

int main(int argc, char *argv[])
{
    //
    // Check for correct arguments
    //
    if (argc < 2)
    {
        std::cout << "Please specify the config file." << std::endl;
        std::cout << "Usage: ./MDMTraceExample <config-file>" << std::endl;
        return 0;
    }

    // create MDMTrace
    MDMTrace* mdmTrace = new MDMTrace();

    //
    // parse json config file
    //
    // create and read json config file
    Json::Value config;
    std::string configFileName = argv[1];
    std::ifstream configStream(configFileName.c_str());
    configStream >> config;
    configStream.close();
    // Read from json config file
    bool usingProbe;
    double mdmAngle;
    double mdmDipoleField;
    double mdmDipoleProbe;
    double mdmMultipoleProbe;
    double scatteredMass;
    double scatteredCharge;
    double scatteredEnergy;
    std::vector<double> scatteredAngles;
    for (Json::Value::iterator it = config.begin(); it != config.end(); it++)
    {
        if (false)
        { // DUMMY - AESTHETICS ONLY
        }
        else if (it.key().asString() == "usingProbe")
        {
            usingProbe = it->asBool();
        }
        else if (it.key().asString() == "mdmAngle")
        { // MDM ANGLE [deg]
            mdmAngle = it->asDouble();
        }
        else if (it.key().asString() == "mdmDipoleField")
        { // MDM DIPOLE FIELD [kG cm]
            mdmDipoleField = it->asDouble();
        }
        else if (it.key().asString() == "mdmDipoleProbe")
        { // MDM DIPOLE PROBE
            mdmDipoleProbe = it->asDouble();
        }
        else if (it.key().asString() == "mdmMultipoleProbe")
        { // MDM MULTIPOLE PROBE
            mdmMultipoleProbe = it->asDouble();
        }
        else if (it.key().asString() == "scatteredMass")
        { // ION MASS [AMU]
            scatteredMass = it->asDouble();
        }
        else if (it.key().asString() == "scatteredCharge")
        { // ION CHARGE STATE [e+ charge]
            scatteredCharge = it->asDouble();
        }
        else if (it.key().asString() == "scatteredEnergy")
        { // ION ENERGY [MeV]
            scatteredEnergy = it->asDouble();
        }
        else if (it.key().asString() == "scatteredAngles")
        { // ION ANGLES [deg]
            for (unsigned int i = 0; i < it->size(); i++)
            {
                scatteredAngles.push_back((*it)[i].asDouble());
            }
        }
    }

    //
    // set MDMTrace parameters
    //
    if(usingProbe) std::cout<<"CONFIRM: Use probe value"<<std::endl;
    if(!usingProbe) std::cout<<"CONFIRM: Not use probe value, use dipole field and 0.71 scale"<<std::endl;
    mdmTrace->SetMDMAngle(mdmAngle);
    std::cout<<"CONFIRM: MDM angle is set to "<<mdmTrace->GetMDMAngle()<<"deg"<<std::endl;
    if(!usingProbe)
    {
        mdmTrace->SetMDMDipoleField(mdmDipoleField);
    }
    else
    {
        mdmTrace->SetMDMProbe(mdmDipoleProbe, mdmMultipoleProbe);
    }
    mdmTrace->SetScatteredMass(scatteredMass);
    std::cout<<"CONFIRM: Scattered mass is set to "<<mdmTrace->GetScatteredMass()<<std::endl;
    mdmTrace->SetScatteredCharge(scatteredCharge);
    std::cout<<"CONFIRM: Scattered charge is set to "<<mdmTrace->GetScatteredCharge()<<std::endl;
    mdmTrace->SetScatteredEnergy(scatteredEnergy);
    std::cout<<"CONFIRM: Scattered energy is set to "<<mdmTrace->GetScatteredEnergy()<<"MeV"<<std::endl;

    for (auto angle : scatteredAngles) 
    {
		mdmTrace->SetScatteredAngle(angle);
		mdmTrace->SendRay();
    double x1,y1,angX1,angY1;
    mdmTrace->GetPositionAngleFirstWire(x1,y1,angX1,angY1);
		//
    printf("Scattered Angle: %.2fdeg  X1: %.2fcm  Y1: %.2fcm  AngX1: %.2fdeg  AngY1: %.2fdeg\n", mdmTrace->GetScatteredAngle(), x1,y1,angX1,angY1);
	}
    return 0;
}