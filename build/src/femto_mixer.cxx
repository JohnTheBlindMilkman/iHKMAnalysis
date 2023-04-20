#include "MixingManager.h"

/*
    _____ _____ __  __ _____ ___    __  __ _____  _______ ____  
    |  ___| ____|  \/  |_   _/ _ \  |  \/  |_ _\ \/ / ____|  _ \ 
    | |_  |  _| | |\/| | | || | | | | |\/| || | \  /|  _| | |_) |
    |  _| | |___| |  | | | || |_| | | |  | || | /  \| |___|  _ < 
    |_|   |_____|_|  |_| |_| \___/  |_|  |_|___/_/\_\_____|_| \_\
                             _____.-._____
                            '-------------'
                            |    (o)(0)   |
                            |  o.(.--).o  |
                            \  O` ) : `o  /
                             | o.( _).O  |
                              \O' `- 'o /
                               |  >|<  |
                               \_______/
                            .'==========='.
                           / o o o o o o o \ 
                          '-----------------'

    Author:         Jędrzej Kołaś
    Date:           09-03-2023
    Affiliation:    Warsaw University of Technology
*/

int main(int argc, char *argv[])
{
    //IOManager::ReadConsole()
    //IOManager::ReadIniFile()
    //IOManager::ReadkTFile()
    //IOManager::InitHistos()
    //MixingManager::Init()
    //IOManager::ReadROOTFiles()

    //MixingManager::Start()

    //IOManager::WriteHistos()
    //MixnigManager::Finish()
    //IOManager::Finish()

    MixingManager MM(argc,argv);
    MM.Init();
    MM.Start();

    return 0;
}