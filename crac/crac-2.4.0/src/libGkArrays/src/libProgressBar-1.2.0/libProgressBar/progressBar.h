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
 * $Id: progressBar.h,v 1.4 2013/11/29 16:46:35 doccy Exp $
 *
 * ---------------------------------------------
 *
 * $Log: progressBar.h,v $
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

#ifndef __PROGRESSBAR_H
#define __PROGRESSBAR_H

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <config.h>

using namespace std;

/**
 * @Namespace DoccY
 * For my personal convenience, I have decided to put all my C++ contributions
 * (at this time, it is the only one...) in this namespace.
 */
namespace DoccY {

  extern "C" {
    /* One can use this function to test the library availability */
    char *libProgressBarVersion();
  }

  /**
   * Class ProgressBar.
   * @Author DoccY <alban.mancheron@inria.fr>
   * @Copyright © 2007 --
   *            Université de Antilles et de la Guyane
   *            Institut National de Recherche en Informatique et en Automatique
   * @brief A C++ class to manage... progress bars.
   */
  class ProgressBar {

  private:

    /**
     * Symbol used for show the progression (default: '=').
     */
    const char symbol;

    /**
     * Starting delimiter (default: '[').
     */
    const char delim1;

    /**
     * Ending delimiter (default: ']').
     */
    const char delim2;

    /**
     * ProgressBar Title.
     */
    const string title;

    /**
     * Current value of the process.
     */
    unsigned int val;

    /**
     * Max value (at the end of the process).
     */
    const unsigned int max_val;

    /**
     * Width of the progress bar.
     * This parameter must be large enough to write the title
     * and the optional time and percents (@see display).
     */
    const unsigned int width;

    /**
     * Output stream to display the progress bar.
     */
    ostream &output;

    /**
     * Starting time of the process.
     */
    clock_t start;

    /**
     * Optional informations to display:
     * if (display & 1), then show percents;
     * if (display & 2), then show consumed time.
     */
    unsigned char display;

  public:
    /**
     * Constructor (standard)
     * @param t ProgressBar title (@see title).
     * @param m Max value at the end of the process (@see max_val).
     * @param w Width of the Progress Bar (@see width).
     * @param os Output stream to display the progress bar (@see output)
     *           [default: ios::cerr].
     * @param startnewline Tell if the progress bar starts by an new line
     *                     [default: true]. This parameter is usefull when
     *                     managing many progress bars.
     * @param start Set starting time to this value (in seconds) [default: 0].
     */
    ProgressBar(const string &t, const unsigned int m, const unsigned int w, ostream &os = cerr, bool startnewline = true, clock_t start = 0);

    /**
     * Constructor (custom)
     * @param s Symbol used for show the progression (@see symbol)
     *          [default: '='].
     * @param d1 Starting delimiter (@see delim1) [default: '['].
     * @param d2 Ending delimiter (@see delim2) [default: '['].
     * @param t ProgressBar title (@see title) [default: "ProgressBar"].
     * @param v Initial value at the begining of the process (@see val)
     *          [default: 0].
     * @param m Max value at the end of the process (@see max_val)
     *          [default: 100].
     * @param w Width of the Progress Bar (@see width) [default: 80].
     * @param os Output stream to display the progress bar (@see output)
     *           [default: ios::cerr].
     * @param startnewline Tell if the progress bar starts by an new line
     *                     [default: true]. This parameter is usefull when
     *                     managing many progress bars.
     * @param start Set starting time to this value (in seconds) [default: 0].
     */
    ProgressBar(const char s = '=', const char d1 = '[', const char d2 = ']',
		const string &t = "ProgressBar",
		unsigned int v = 0, const unsigned int m = 100, const unsigned int w = 80,
		ostream &os = cerr, bool startnewline = true, clock_t start = 0);

    /**
     * Destructor
     */
    ~ProgressBar();

    /**
     * Enable percentage display (@see display).
     */
    void ShowPercent();

    /**
     * Disable percentage display (@see display).
     */
    void HidePercent();

    /**
     * Enable consumed time display (@see display).
     */
    void ShowTime();

    /**
     * Disable consumed time display (@see display).
     */
    void HideTime();

    /**
     * Get consumed time in seconds (@see start).
     */
    double GetTime() const;

    /**
     * Add sec seconds to expanded time (@see start).
     */
    void AddTime(clock_t sec);

    /**
     * Reset the progress bar (@see val, @see start).
     */
    void Reset();

    /**
     * Set the current value of the process.
     * @param v The current value (@see val).
     * @param update Refresh the progress bar display.
     */
    void SetVal(unsigned int v, bool update = false);

    /**
     * Increment by one the current value of the process
     * (@see val).
     * @param update Refresh the progress bar display.
     */
    void Step(bool update = true);

    /**
     * Refresh the progress bar display.
     * @param goback If true, rewind the cursor at the
     *               beginning of the progress bar.
     */
    void update(bool goback = true) const;

    /**
     * Clear the progress bar from the display. This
     * doesn't reset either the timer, or the current
     * value!!!
     */
    void clear() const;

    /**
     * Fill the progress bar with some symbol.
     * @param c Symbol used to fill the progress bar
     *                 [default: ' '].
     */
    void fill(const char c = ' ') const;
  }; /* Fin ProgressBar */

} /* Fin espace de nom DoccY */

#endif /*  __PROGRESSBAR_H */
// Local Variables:
// mode:c++
// End:
