#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QPushButton>
#include <QTextEdit>
#include <QFileDialog>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QString>
#include <QTextStream>
#include <QRegularExpression>
#include <QMessageBox>

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void addSequence();
    void searchKmer();
    void predictFunction();
    void showGCContent();

private:
    QTextEdit *resultWindow;
    QString currentSequence;
    double calculateGC(const QString &sequence);
};

#endif // MAINWINDOW_H

