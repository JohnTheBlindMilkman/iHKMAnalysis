/********************************************************************************
 *                                                                              *
 *             THERMINATOR 2: THERMal heavy-IoN generATOR 2                     *
 *                                                                              *
 * Version:                                                                     *
 *      Release, 2.0.3, 1 February 2011                                         *
 *                                                                              *
 * Authors:                                                                     *
 *      Mikolaj Chojnacki   (Mikolaj.Chojnacki@ifj.edu.pl)                      *
 *      Adam Kisiel         (kisiel@if.pw.edu.pl)                               *
 *      Wojciech Broniowski (Wojciech.Broniowski@ifj.edu.pl)                    *
 *      Wojciech Florkowski (Wojciech.Florkowski@ifj.edu.pl)                    *
 *                                                                              *
 * Project homepage:                                                            *
 *      http://therminator2.ifj.edu.pl/                                         *
 *                                                                              *
 * For the detailed description of the program and further references           *
 * to the description of the model please refer to                              *
 * http://arxiv.org/abs/1102.0273                                               *
 *                                                                              *
 * This code can be freely used and redistributed. However if you decide to     *
 * make modifications to the code, please, inform the authors.                  *
 * Any publication of results obtained using this code must include the         *
 * reference to arXiv:1102.0273 and the published version of it, when           *
 * available.                                                                   *
 *                                                                              *
 ********************************************************************************/

#ifndef _TH2_GLOBAL_H_
  #define _TH2_GLOBAL_H_

#include <iostream>
#include <stdlib.h>
#include <cstdio>
#include <TString.h>

namespace THGlobal
{
    // Global static constants
    static const double kTwoPi3 = 248.0502134423985614038105205368;
    static const double kTwoPi2 = 39.4784176043574344753379639995;
    static const double kHbarC = 0.1973269631;
    static const float GeVtoFm = 0.197327;

    static const int kFileNameMaxChar = 2000;
    static const int NoParams = 7;
    static const int projWidth[] = {0,4,10};

    static const TString sParNames[NoParams] = {"LambdaInv","Rinv","Norm","LambdaOSL","Rout","Rside","Rlong"};
    static const TString sProjNames[] = {"out","side","long","inv"};

    static const TString sTreeName = "treefin";

    enum ReferenceFrame
    {
        LCMS = 1,
        PCMS = 2,
        PRF = 3
    };

    enum INIParam
    {
        EisPrimordial,
        EuseCoulomb,
        EwriteHisto,
        EpairType,
        EeventToMix,
        EeventFiles,
        EtimeCut,
        EktMin,
        EktMax
    };

    // Define version of THERMINATOR 2
    static const TString therminatorVersion = "2.1.1";

    // Define compilation specific variables
    static const TString cxxVersion = Form("g++(%d.%d.%d)", __GNUC__ ,__GNUC_MINOR__ ,__GNUC_PATCHLEVEL__);
    static const TString rootVersion = ". / ";

    // Define global errors
    enum Error
    {
        generalNoError,
        generalFileNotFound, 
        generalMidelUnknown,
        generalunsuportedValue,
        generalUnknownRefFrame,
        configParameterNotFound,
        libraryEmpty,
        libraryTagNotFound,
        libraryTagAtribNotFound,
        libraryVector3DNotFound,
        femtoUnknownPairType,
        femtoWrongkTOrder,
        femtoGraphNotFound
    };

    // Define particle structure
    struct ParticleCoor
    {
        float px,py,pz,e,x,y,z,t;
        int id,mid;
    };

    // Define DEBUG information
    #define PRINT_MESSAGE(_mes) std::cout << _mes << std::endl;
    #if _DEBUG_LEVEL_==0
        #define PRINT_DEBUG_3(_mes)
        #define PRINT_DEBUG_2(_mes)
        #define PRINT_DEBUG_1(_mes)
    #elif _DEBUG_LEVEL_==1
        #define PRINT_DEBUG_3(_mes)
        #define PRINT_DEBUG_2(_mes)
        #define PRINT_DEBUG_1(_mes) std::cerr << _mes << std::endl;
    #elif _DEBUG_LEVEL_==2
        #define PRINT_DEBUG_3(_mes)
        #define PRINT_DEBUG_2(_mes) std::cerr << _mes << std::endl;
        #define PRINT_DEBUG_1(_mes) std::cerr << _mes << std::endl;
    #elif _DEBUG_LEVEL_==3
        #define PRINT_DEBUG_3(_mes) std::cerr << "DL3:" << _mes << std::endl;
        #define PRINT_DEBUG_2(_mes) std::cerr << "DL2:" << _mes << std::endl;
        #define PRINT_DEBUG_1(_mes) std::cerr << "DL1:" << _mes << std::endl;
    #endif
}

#endif /*_TH2_GLOBAL_H_*/


/*! @file THGlobal.h
 * @brief <c><b>THERMINATOR 2</b></c> global definitions.
 */
/*! @var static const double kTwoPi3
 * @brief @f$ \left(2 \pi\right)^3 @f$
 *
 * @var static const double kTwoPi2
 * @brief @f$ \left(2 \pi\right)^2 @f$
 *
 * @var static const double kHbarC
 * @brief @f$ \hbar c @f$ = 0.1973269631 [GeV fm]
 */
/*! @def _THERMINATOR2_VERSION_
 * @brief Version number of <c><b>THERMINATOR 2</b></c>
 */
/*! @def _CXX_VER_
 * @brief Version of the C++ compiler.
 *
 * @def _ROOT_VER_
 * @brief Version of ROOT.
 */
/*! @def _ERROR_GENERAL_FILE_NOT_FOUND_
 * @brief Integer returned if file was not found.
 *
 * @def _ERROR_GENERAL_MODEL_UNKNOWN_
 * @brief Integer returned if Model name is unknown.
 *
 * @def _ERROR_CONFIG_PARAMETER_NOT_FOUND_
 * @brief Integer returned if Configurator will not find Keyword. 
 *
 * @def _ERROR_LIBRARY_EMPTY_
 * @brief Integer returned by Hypersurface_Library if library is empty.
 *
 * @def _ERROR_LIBRARY_TAG_NOT_FOUND_
 * @brief Integer returned by Hypersurface_Library if XML tag is not found.
 *
 * @def _ERROR_LIBRARY_TAG_ATTRIB_NOT_FOUND_
 * @brief Integer returned by Hypersurface_Library if XML tag attribute is not found.
 *
 * @def _ERROR_LIBRARY_VECTOR3D_NOT_FOUND_
 * @brief Integer returned by Hypersurface_Library if Vector3D information is not found.
 *
 * @def _ERROR_FEMTO_UNKNOWN_PAIRTYPE_
 * @brief Integer returned by therm_femto if pair type is not known.
 */
/*! @def PRINT_MESSAGE(_mes)
 * @brief Part of the debugging system. Sends <em>_mes</em> to std::cout. 
 *
 * @def PRINT_DEBUG_1(_mes)
 * @brief Part of the debugging system. Sends <em>_mes</em> to std::cerr if the compilation variable <em>_DEBUG_LEVEL_ >= 1</em>
 *
 * @def PRINT_DEBUG_2(_mes)
 * @brief Part of the debugging system. Sends <em>_mes</em> to std::cerr if the compilation variable <em>_DEBUG_LEVEL_ >= 2</em>
 *
 * @def PRINT_DEBUG_3(_mes)
 * @brief Part of the debugging system. Sends <em>_mes</em> to std::cerr if the compilation variable <em>_DEBUG_LEVEL_ = 3</em>
 */
