// Author(s): Ruben Vink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "graphicsinfodialog.h"


GraphicsInfoDialog::GraphicsInfoDialog(QWidget *parent):
  QDialog(parent)
{
  m_ui.setupUi(this);
  lbl_info = m_ui.lbl_info;
}