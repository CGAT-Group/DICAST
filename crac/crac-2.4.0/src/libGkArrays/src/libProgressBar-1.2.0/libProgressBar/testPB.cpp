/******************************************************************************
*                                                                             *
*     Copyright © 2007      -- Université de Antilles et de la Guyane         *
*     Copyright © 2007      -- Institut National de Recherche en Informatique *
*                              et en Automatique                              *
*                                                                             *
* Author: DoccY <alban.mancheron@inria.fr>                                    *
*                                                                             *
* This File initially comes from StatiSTARS.                                  *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Library General Public License as published by   *
* the Free Software Foundation; either version 2 of the License, or (at your  *
* option) any later version.                                                  *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public       *
* License for more details.                                                   *
*                                                                             *
* You should have received a copy of the GNU Library General Public License   *
* along with this library; if not, write:                                     *
*                                                                             *
*                    the Free Software Foundation, Inc.,                      *
*                    59 Temple Place - Suite 330,                             *
*                    Boston, MA  02111-1307, USA.                             *
*                                                                             *
******************************************************************************/

/*
 * =============================================
 *
 * $Id: testPB.cpp,v 1.3 2013/11/29 16:46:35 doccy Exp $
 *
 * ---------------------------------------------
 *
 * $Log: testPB.cpp,v $
 * Revision 1.3  2013/11/29 16:46:35  doccy
 * Update:
 * - Time display enhancement.
 * - Allows systems without tex distribution to install the library (doc is just not built).
 * - Add a method/constructor option to add/set an initial  a value to the timer.
 * - Add a 'C' function to get version number.
 *
 * Revision 1.2  2007/09/27 14:40:50  mancheron
 * Ajout de la méthode GetTime() qui permet de récupérer le temps
 * écoulé (en s.) depuis le début du processus en progrès.
 *
 * =============================================
 *
 */

#include "progressBar.h"

#define MAX_VAL 1000

using namespace DoccY;

void testPB(ProgressBar &pb, double tps, bool fake = false) {

  pb.Reset();
  double ftps = tps / MAX_VAL;
  if (fake) {
    tps = 30;
  }
  tps /= MAX_VAL;
  for(unsigned int i = 0; i < MAX_VAL; i++) {
    clock_t timer = clock();
    pb.Step();
    if (fake) {
      pb.AddTime(ftps);
    }
    while ((clock() - timer)/double(CLOCKS_PER_SEC) < tps);
  } /* Fin Pour */
  if (fake) {
    pb.Reset();
    pb.AddTime(ftps * MAX_VAL);
    pb.SetVal(MAX_VAL);
  }
  pb.update(false);
  cout << endl;

} /* Fin testPB */

int main(int argc, char** argv) {

  ProgressBar PB("ProgressBar Test Program", MAX_VAL, 80, cout, true);
  for (int i = 1; i <= 8; i++) {
    cout << setfill(' ') << setw(10) << i*10;
  }
  cout << endl;
  for (int i = 0; i < 8; i++) {
    cout << setfill('-') << setw(10) << "+";
  }
  cout << endl;

  cout << "Test with percent and without time [default] (10s):" << endl;
  testPB(PB, 2.5);

  cout << "Test with percent and time (10s):" << endl;
  PB.ShowTime();
  testPB(PB, 2.5);

  cout << "Test without percent and with time (10s):" << endl;
  PB.HidePercent();
  testPB(PB, 2.5);

  cout << "Test without percent and without time (10s)" << endl;
  PB.HideTime();
  testPB(PB, 2.5);

  cout << "Test without percent and with time (2 min)" << endl;
  PB.ShowPercent();
  PB.ShowTime();
  testPB(PB, 120.0);

  cout << "Test without percent and with time (1h30 [displayed in 30s])" << endl;
  testPB(PB, (1*60.0+30.0)*60.0, true);
  cout << "Test without percent and with time (15h15 [displayed in 30s])" << endl;
  testPB(PB, (15*60.0+15.0)*60.0, true);
  cout << "Test without percent and with time (22d22 [displayed 30s])" << endl;
  testPB(PB, (22.0*24.0+22.0)*3600.0, true);

  return 0;
}
