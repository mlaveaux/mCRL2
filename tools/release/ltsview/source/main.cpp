// Author(s): Bas Ploeger
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//


#include "mcrl2/utilities/input_tool.h"
#include "mcrl2/gui/qt_tool.h"
#include "mainwindow.h"

using namespace mcrl2;
using namespace mcrl2::utilities;

typedef gui::qt::qt_tool<tools::input_tool> ltsview_base;
class ltsview_tool : public ltsview_base
{
  public:
    ltsview_tool():
      ltsview_base("LTSView",
        (QStringList() << "Bas Ploeger" << "Carst Tankink" << "Ruud Koolen").join(", ").toStdString(),
        "3D interactive visualization of a labelled transition system", // what-is
        "Start the LTSView application. If INFILE is supplied then the " // command-line description
        "LTS in INFILE is loaded into the application.\n"
        "\n"
        "The input format is determined by the contents of INFILE. If that fails, "
        "an attempt is made to force the input format based on the file extension. "
        "The supported formats with file extensions are:\n"
        "  Aldebaran format (CADP; *.aut);\n"
        "  GraphViz format (*.dot);\n"
        "  Finite State Machine format (*.fsm);\n"
        "  mCRL SVC format (*.svc);\n"
        "  mCRL2 format (*.lts).",
        "Tool for interactive visualization of state transition systems.\n" // gui-specific description
        "\n"
        "LTSView is based on visualization techniques by Frank van Ham and Jack van Wijk.\n"
        "See: F. van Ham, H. van de Wetering and J.J. van Wijk,\n"
        "\"Visualization of State Transition Graphs\". "
        "Proceedings of the IEEE Symposium on Information Visualization 2001. IEEE CS Press, pp. 59-66, 2001.\n"
        "\n"
        "The default colour scheme for state marking was obtained from http://www.colorbrewer.org",
        "http://www.mcrl2.org/web/user_manual/tools/release/ltsview.html")
    {}

    bool run()
    {
      qRegisterMetaType<LTS *>("LTS *");

      QThread atermThread;
      atermThread.start();

      MainWindow window(&atermThread);

      if (!m_input_filename.empty())
      {
        QMetaObject::invokeMethod(&window, "open", Qt::QueuedConnection,
                           Q_ARG(QString, QString::fromStdString(m_input_filename)));
      }

      return show_main_window(window);
    }
};

int main(int argc, char *argv[])
{
  return ltsview_tool().execute(argc, argv);
}
