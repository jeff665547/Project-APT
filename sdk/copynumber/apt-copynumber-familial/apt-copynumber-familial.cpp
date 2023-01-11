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
 * @file   apt-copynumber-familial.cpp
 *
 * @brief  Initial version of program that uses copynumber architecture to
 *         analyze trios and duos or cychp files.
 */

//
#include "copynumber/CNFamilialAnalysisMethodFactory.h"
#include "copynumber/CNFamilialEngine.h"
//
#include "chipstream/EngineUtil.h"
#include "util/AptVersionInfo.h"
#include "util/Fs.h"
#include "util/LogStream.h"
#include "util/PgOptions.h"
#include "util/Util.h"

//
#ifdef _DEBUG
#ifndef _WIN64
#ifdef _USE_VLD
#include "../external/vld/vld.h"
#endif
#endif
#endif

#ifndef WIN32
#include <unistd.h>
#endif /* WIN32 */


/** Everybody's favorite function... */
int main(int argc,const char* argv[]) {
  ofstream logOut;
  LogStream log;
  string logName;

  try {

    CNFamilialEngine engine;

    const string version = AptVersionInfo::version();
    const string cvsId = AptVersionInfo::cvsId();
    const string versionToReport = AptVersionInfo::versionToReport();
    const string execGuid = affxutil::Guid::GenerateNewGuid();

    /* Parse options. */
    engine.setUsage(
         "apt-copynumber-familial - program for analyzing trios and duos of cychp files\n"
         "from Affymetrix cytogenetics microarrays.\n"
         "\n"
         "usage:\n"
         "./apt-copynumber-familial \\\n"
           "     --index-cychp-file Child.Cytogenetics.cychp  \\\n"
           "     --mother-cychp-file Mother.Cytogenetics.cychp  \\\n"
           "     --father-cychp-file Father.Cytogenetics.cychp \\\n"
           "     --allele-frequency-file AlleleFrequencies_Reduced.txt \\\n"
           "     --analysis paternity.role-validity-threshold=1000.0 \\\n"
           "     --out-dir outputDir \\\n"
           "     --familial-file Trio_1_1_3_2.fam");

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

    engine.openStandardLog("apt-copynumber-familial.log",logOut,log);

    engine.run();
  }
  catch(const Except& ex) {
      Verbose::out(1, ex.what());
      return 1;
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
