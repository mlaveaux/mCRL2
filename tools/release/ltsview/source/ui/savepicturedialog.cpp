// Author(s): Bas Ploeger, Carst Tankink, Ruud Koolen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "savepicturedialog.h"
#include <QImageWriter>

SavePictureDialog::SavePictureDialog(QWidget *parent, LtsCanvas *canvas, QString filename):
  QDialog(parent),
  m_canvas(canvas),
  m_filename(filename),
  m_inChange(false)
{
  m_ui.setupUi(this);

  m_width = canvas->viewWidth();
  m_height = canvas->viewHeight();
  m_ui.width->setValue(m_width);
  m_ui.height->setValue(m_height);

  connect(m_ui.width, SIGNAL(valueChanged(int)), this, SLOT(widthChanged(int)));
  connect(m_ui.height, SIGNAL(valueChanged(int)), this, SLOT(heightChanged(int)));
  connect(m_ui.buttonBox, SIGNAL(accepted()), this, SLOT(save()));
  connect(m_ui.buttonBox, SIGNAL(rejected()), this, SLOT(reject()));
}

void SavePictureDialog::widthChanged(int value)
{
  if (m_inChange)
  {
    return;
  }

  if (m_ui.maintainAspectRatio->isChecked())
  {
    m_inChange = true;
    m_ui.height->setValue((int)(value * m_height / m_width));
    m_inChange = false;
  }
}

void SavePictureDialog::heightChanged(int value)
{
  if (m_inChange)
  {
    return;
  }

  if (m_ui.maintainAspectRatio->isChecked())
  {
    m_inChange = true;
    m_ui.height->setValue((int)(value * m_width / m_height));
    m_inChange = false;
  }
}

void SavePictureDialog::save()
{
  int width = m_ui.width->value();
  int height = m_ui.height->value();

  QImage image = m_canvas->renderImage(width, height);

  emit statusMessage("Saving image...");
  QImageWriter writer(m_filename);
  if (writer.write(image))
  {
    emit statusMessage("Done");
  }
  else
  {
    emit statusMessage("Saving image failed.");
  }

  accept();
}
