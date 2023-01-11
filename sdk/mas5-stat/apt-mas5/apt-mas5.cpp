////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify 
// it under the terms of the GNU General Public License (version 2) as 
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License 
// along with this program;if not, write to the 
// 
// Free Software Foundation, Inc., 
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////

/**
 * @file   apt-mas5.cpp
 * @author Alan Williams
 * @date   Wed Apr 15 15:29:23 PDT 2009
 *
 * @brief  wrapper around engine to run MAS5 algos
 */

//
#include "mas5-stat/apt-mas5/MAS5Engine.h"
//
#include "chipstream/EngineUtil.h"
#include "util/AptVersionInfo.h"
#include "util/BaseEngine.h"
#include "util/Fs.h"
#include "util/LogStream.h"
#include "util/PgOptions.h"
#include "util/Util.h"

//
#ifndef WIN32
#include <unistd.h>
#endif


/** Everybody's favorite function... */
int main(int argc,const char* argv[]) {
  ofstream logOut;
  LogStream log;
  string logName;
    
  try {
    MAS5Engine engine;
  
    const string version = AptVersionInfo::version();
    const string cvsId = AptVersionInfo::cvsId();
    const string versionToReport = AptVersionInfo::versionToReport();
    const string execGuid = affxutil::Guid::GenerateNewGuid();
    
    /* Parse options. */
    engine.setUsage(
                    "apt-mas5 - a program to run the MAS5 algoriths\n"
                    "\n"
                    "usage:\n"
                    "    apt-mas5 .... \n"
                    "          *.cel\n"
                    );
    engine.parseArgv(argv);
    Verbose::setLevel(engine.getOptInt("verbose"));
    const string progName = Fs::basename(engine.getProgName());
    engine.setOpt("command-line",engine.commandLine());
    engine.setOpt("program-name",progName);
    engine.setOpt("program-company","Affymetrix");
    engine.setOpt("program-version",version);
    engine.setOpt("program-cvs-id",cvsId);
    engine.setOpt("version-to-report",versionToReport);
    engine.setOpt("exec-guid",execGuid);
    if(argc == 1) { engine.setOpt("help","true"); }
    
    // Check Options. Will print out version/help if requested then exit.
    engine.checkOptions();
    
    /* Set up the logging and message handlers. */
    if(!Fs::isWriteableDir(engine.getOpt("out-dir"))) {
        if(Fs::mkdirPath(engine.getOpt("out-dir"), false) != APT_OK) {
            Err::errAbort("Can't make or write to directory: " + engine.getOpt("out-dir"));
        }
    }
    if(engine.getOpt("log-file") != "")
        logName = engine.getOpt("log-file");
    else
      logName = Fs::join(engine.getOpt("out-dir") , ToStr("apt-mas5.log"));
    Fs::mustOpenToWrite(logOut, logName.c_str());
    log.setStream(&logOut);
    Verbose::pushMsgHandler(&log);
    Verbose::pushProgressHandler(&log);
    Verbose::pushWarnHandler(&log);
    
    engine.run();
  } 
  catch(...) {
      Verbose::out(1,"Unexpected Error: uncaught exception.");
      // Close log files
      logOut.close();
      return 1;
  }
  // Close log files
  logOut.close();

  return 0;
}
