﻿#include <IQWidgets/igQtColorBarWidget.h>
#include <QDebug>
#include <iGameSceneManager.h>
#include "iGameSmartPointer.h"
#include "iGameDrawObject.h"
//#include <iGameModelColorManager.h>
#include <QMouseEvent>
#include <iomanip>
igQtColorBarWidget::igQtColorBarWidget(QWidget* parent) : QWidget(parent)
{
	initDrawStringStyle();

	this->setWindowFlags(Qt::Widget | Qt::WindowStaysOnTopHint);
	this->setWindowOpacity(0.5);
	//this->setWindowFlags(Qt::FramelessWindowHint);
	//this->setAttribute(Qt::WA_TranslucentBackground,true);
	this->limitX1 = 0;
	this->limitY1 = 0;
	this->limitX2 = 10000;
	this->limitY2 = 10000;
	this->limitHeight = 200;
	this->setMaximumHeight(limitHeight);
	this->resize(this->width(), limitHeight);
	this->move(0, this->limitY2 - 400);

}
void igQtColorBarWidget::updateColorBarDrawInfo() {

    auto scene = iGame::SceneManager::Instance()->GetCurrentScene();
    m_ColorMapper = nullptr;
	if (scene) {
        auto model = scene->GetCurrentModel();
        if (model&&model->GetDataObject()) {
			m_ColorMapper =model->GetDataObject()->GetColorMapper();
		}
	}
    if (!m_ColorMapper) { 
		m_ColorMapper = iGame::ScalarsToColors::New();
	}

	// 计算colors
	this->colorBarLength = this->height() - this->fontHeight - 1-10;
	colors.resize(colorBarLength);
	float scale = 1.0 / colorBarLength;
	QColor color;
	float rgb[16];
	for (int i = 0; i < colorBarLength; i++) {
		float val = i * scale;
        m_ColorMapper->MapColor(val, rgb);
		color.setRgbF(rgb[0], rgb[1], rgb[2]);
		colors[i] = (color);
	}
    float min = m_ColorMapper->GetRange()[0];
    float max = m_ColorMapper->GetRange()[1];
	this->strData.resize(6);
	for (int i = 0; i < 6; i++) {
		float x = min + i * (max - min) * 0.2;
		strData[i] = QString::fromStdString(formatNumber(x));
	}

}
std::string igQtColorBarWidget::formatNumber(float x)
{
	bool isNegative = false;
	std::ostringstream strStream;
	int precision = 4;
	if (x < 0) {
		x = -x;
		isNegative = true;
	}
	std::string str = isNegative == true ? "-" : "";
	if (x < 0.01 && x >= 0.0000000001) {
		
		float a = x;
		int num = 0;
		precision = 3;
		while (a < 1) {
			a *= 10;
			num++;
		}
		strStream << std::setprecision(precision) << a;
		str += strStream.str() + "e-" + std::to_string(num);
	}
	else if (x > 1000) {
		float a = x;
		int num = 0;
		while (a > 10) {
			a /= 10.0;
			num++;
		}
		strStream << std::setprecision(precision) << a;
		str += strStream.str() + 'e' + std::to_string(num);
	}
	else {
		if (x < 0.1)precision = 4;
		else if (x < 1)precision = 5;
		else precision = 6;
		strStream << std::setprecision(precision) << x;
		str += strStream.str();
	}
	return str;
    //return "";
}

void igQtColorBarWidget::initDrawStringStyle()
{
	this->font = *(new QFont);
	font.setFamily("Microsoft YaHei");
	font.setPointSize(8);
	QFontMetrics fm(this->font);
	this->fontHeight = fm.boundingRect('1').height();
}
void igQtColorBarWidget::mousePressEvent(QMouseEvent* _event)
{
	this->isPressed = true;
	this->boundColor = Qt::white;
	this->lastPos = _event->pos();
}
void igQtColorBarWidget::mouseMoveEvent(QMouseEvent* _event)
{
	auto point = _event->pos();
	int offsetX = point.x() - lastPos.x();
	int offsetY = point.y() - lastPos.y();
	auto x = this->x() + offsetX;
	auto y = this->y() + offsetY;
	x = qBound(limitX1, x, limitX2 - this->width());
	y = qBound(limitY1, y, limitY2 - this->height());
	this->move(x, y);
}
void igQtColorBarWidget::mouseReleaseEvent(QMouseEvent* _event)
{
	this->isPressed = false;
	this->boundColor = Qt::black;
	update();

}
void igQtColorBarWidget::setLimitPosition(int x1, int y1, int x2, int y2)
{
	this->limitX1 = x1;
	this->limitX2 = x2;
	this->limitY1 = y1;
	this->limitY2 = y2;

	auto x = this->x();
	auto y = this->y();
	x = qBound(limitX1, x, limitX2 - this->width());
	y = qBound(limitY1, y, limitY2 - this->height());
	this->move(x, y);
}

void igQtColorBarWidget::paintEvent(QPaintEvent* event)
{
	updateColorBarDrawInfo();
	QPainter painter(this);
	QVector<QRect>data;
	data.resize(6);
	QFontMetrics fm(this->font);
	int maxTextWidth = -1;
	for (int i = 0; i < 6; i++) {
		QRect rec = fm.boundingRect(strData[i]);
		data[i] = rec;
		maxTextWidth = std::max(maxTextWidth, rec.width());
	}
	int height = this->height();
	int width = this->width();




	//draw colorbar
	int minColorBarWidth = 30;
	int maxColorBarWidth = 150;
	//this->colorBarWidth = std::max(width - maxTextWidth - 20, minColorBarWidth);
	this->colorBarWidth = minColorBarWidth;
	int st = 5+colorBarLength + this->fontHeight / 2;
	for (int i = 0; i < colorBarLength; i++) {
		QRect rect(0, st - i, colorBarWidth, 1);
		painter.fillRect(rect, colors.at(i));
	}
	painter.fillRect(QRect(0, st - colorBarLength, 2, colorBarLength), this->boundColor);
	painter.fillRect(QRect(colorBarWidth - 2, st - colorBarLength, 2, colorBarLength), this->boundColor);
	painter.fillRect(QRect(0, st - colorBarLength, colorBarWidth, 2), this->boundColor);
	painter.fillRect(QRect(0, st - 2, colorBarWidth, 2), this->boundColor);


	// draw datarange
	painter.setFont(this->font);
	int restWidth = 5;
	for (int i = 0; i < data.size(); i++) {
		int textWidth = data[i].width();
		int textHeight = data[i].height();
		float y = st * 1.0 - i * 1.0 * colorBarLength / 5.0;
		QRect textRect(colorBarWidth + 2, (int)y - (textHeight + 1) / 2, textWidth + restWidth, textHeight);
		painter.drawText(textRect, strData[i]);
		if (i == 0 || i == data.size() - 1)continue;
		painter.fillRect(QRect(2, (int)y, colorBarWidth / 4, 2), this->boundColor);
		painter.fillRect(QRect(3 * colorBarWidth / 4, (int)y, colorBarWidth - 3 * colorBarWidth / 4, 2), this->boundColor);

	}
	int minWidth = maxTextWidth + minColorBarWidth + restWidth;
	//int maxWidth = maxTextWidth + maxColorBarWidth + restWidth;
	//this->setMinimumWidth(minWidth);
	//this->setMaximumWidth(maxWidth);
	//this->setMaximumWidth(minWidth);
	this->resize(minWidth, this->limitY2 - 5);
	//std::cout << minWidth << " " << this->width() << '\n';
	Q_EMIT PaintFinished();
}

