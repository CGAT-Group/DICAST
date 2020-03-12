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
 * $Id: progressBar.cpp,v 1.4 2013/11/29 16:46:35 doccy Exp $
 *
 * ---------------------------------------------
 *
 * $Log: progressBar.cpp,v $
 * Revision 1.4  2013/11/29 16:46:35  doccy
 * Update:
 * - Time display enhancement.
 * - Allows systems without tex distribution to install the library (doc is just not built).
 * - Add a method/constructor option to add/set an initial  a value to the timer.
 * - Add a 'C' function to get version number.
 *
 * Revision 1.3  2007/09/27 14:40:50  mancheron
 * Ajout de la méthode GetTime() qui permet de récupérer le temps
 * écoulé (en s.) depuis le début du processus en progrès.
 *
 * Revision 1.2  2007-09-09 21:16:29  mancheron
 * Modulation de l'affichage (avec ou sans temps/pourcentage).
 *
 * Revision 1.1  2007-08-15 16:36:23  mancheron
 * *** empty log message ***
 *
 * =============================================
 *
 */

#include "progressBar.h"

using namespace DoccY;

extern "C" {
  char *libProgressBarVersion() {
    return (char*)PACKAGE_VERSION;
  }
}

void ProgressBar::ShowPercent() {
  display |= 1;
} /* Fin ShowPercent */

void ProgressBar::HidePercent() {
  display &= (unsigned char) -2;
} /* Fin HidePercent */

void ProgressBar::ShowTime() {
  display |= 2;
} /* Fin ShowTime */

void ProgressBar::HideTime() {
  display &= (unsigned char) -3;
} /* Fin HideTime */

double ProgressBar::GetTime() const {
  return (clock() - start)/double(CLOCKS_PER_SEC);
} /* Fin GetTime */

void ProgressBar::AddTime(clock_t sec) {
  start -= sec * CLOCKS_PER_SEC;
}

void ProgressBar::Reset() {
  val = 0;
  start = clock();
} /* Fin Reset */

void ProgressBar::SetVal(unsigned int v, bool update) {
  val = min(v, max_val);
  if (update) {
    this->update();
  } /* Fin Si */
} /* Fin SetVal */

void ProgressBar::Step(bool update) {
  val = min(val+1, max_val);
  if (update) {
    this->update();
  } /* Fin Si */
} /* Fin Step */

ProgressBar::ProgressBar(const string &t, const unsigned int m, const unsigned int w, ostream &os, bool startnewline, clock_t start):
  symbol('='), delim1('['), delim2(']'), title(t), val(0), max_val(m),
  width(w), output(os), start(clock() - start*CLOCKS_PER_SEC), display(1) {
  if (startnewline) {
    output << endl;
  } /* Fin Si */
} /* Fin ProgressBar */

ProgressBar::ProgressBar(const char s, const char d1, const char d2, const string &t,
			      unsigned int v, const unsigned int m, const unsigned int w, 
			      ostream &os, bool startnewline, clock_t start):
  symbol(s), delim1(d1), delim2(d2), title(t), val(min(v, m)), max_val(m),
  width(w), output(os), start(clock() - start*CLOCKS_PER_SEC), display(1) {
  if (startnewline) {
    output << endl;
  } /* Fin Si */
} /* Fin ProgressBar */

ProgressBar::~ProgressBar() {
  clear();
} /* Fin ~ProgressBar */

void ProgressBar::update(bool goback) const {
  unsigned int lg;
  unsigned int v;
  stringstream tps;
  if (display & 1) {
    tps << setw(4) << setfill(' ') << val * 100 / max_val << " %";
  } /* Fin Si */
  if (display & 2) {
    double temps = (clock() - start)/double(CLOCKS_PER_SEC);
    if (display & 1) {
      tps << " (";
    } else {
      tps << " ";
    }
    if (temps < 60.0) {
      tps.precision((temps < 10.0)?2:1);
      tps.setf(stringstream::fixed,stringstream::floatfield);
      tps << temps << "s";
    } else {
      tps.precision(0);
      if (temps < 3600.0) {
	tps << int(temps)/60 << "m"
	    << setfill('0') << setw(2) << int(temps)%60;
      } else {
	if (temps < 86400.0) {
	  tps << int(temps)/3600 << "h"
	      << setfill('0') << setw(2) << int(temps)%3600/60;
	} else {
	  tps << int(temps)/86400 << "d"
	      << setfill('0') << setw(2) << int(temps)%86400/3600;
	} /* Fin Si */
      } /* Fin Si */
    } /* Fin Si */
    if (display & 1) {
      tps << ")";
    }
  } /* Fin Si */
  lg = width - title.size() - 4 - tps.str().length();
  v = val * lg / max_val;
  output << title << " " << delim1;
  if (v) {
    output << setw(v) << setfill(symbol) << symbol;
  } /* Fin Si */
  output << setw(lg - v + 1) << setfill(' ') << delim2;
  output << tps.str() << flush;
  if (goback) {
    fill('\b');
  } /* Fin Si */
} /* Fin update */

void ProgressBar::clear() const {
  fill(' ');
  fill('\b');
} /* Fin clear */

void ProgressBar::fill(const char c) const {
  output << setw(width - 1) << setfill(c) << c << flush;
} /* Fin fill */
