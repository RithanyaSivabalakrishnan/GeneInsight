#include "mainwindow.h"
#include "dna_analyzer.h"
#include <QFile>
#include <QDebug>
#include <QInputDialog>
#include <QLineEdit>
#include <QTextStream>
#include <QMessageBox>
#include <QFileDialog>
#include <QLabel>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    // Window title
    setWindowTitle("DNA Analyzer 🧬");

    QWidget *centralWidget = new QWidget(this);
    QVBoxLayout *mainLayout = new QVBoxLayout();

    // --- Top Menu Buttons ---
    QHBoxLayout *menuLayout = new QHBoxLayout();

    QPushButton *addBtn = new QPushButton("🧬 Add");
    QPushButton *searchBtn = new QPushButton("🔍 Search");
    QPushButton *predictBtn = new QPushButton("🧫 Predict");
    QPushButton *gcBtn = new QPushButton("📊 GC Content");

    QString btnStyle =
        "background-color: #1e90ff; "
        "color: white; font-weight: bold; "
        "padding: 8px 12px; border-radius: 8px; "
        "font-size: 14px;";

    addBtn->setStyleSheet(btnStyle);
    searchBtn->setStyleSheet(btnStyle);
    predictBtn->setStyleSheet(btnStyle);
    gcBtn->setStyleSheet(btnStyle);

    menuLayout->addWidget(addBtn);
    menuLayout->addWidget(searchBtn);
    menuLayout->addWidget(predictBtn);
    menuLayout->addWidget(gcBtn);

    // --- Result Display Window ---
    resultWindow = new QTextEdit();
    resultWindow->setReadOnly(true);
    resultWindow->setMinimumHeight(450);
    resultWindow->setStyleSheet(
        "border: 2px solid #1e90ff; "
        "border-radius: 10px; "
        "background-color: white; "
        "color: black; "
        "font-family: Consolas; font-size: 14px; "
        "padding: 10px;"
        );

    mainLayout->addLayout(menuLayout);
    mainLayout->addWidget(new QLabel("----------------------------------------------------------"));
    mainLayout->addWidget(resultWindow);
    mainLayout->addWidget(new QLabel("----------------------------------------------------------"));

    centralWidget->setLayout(mainLayout);
    setCentralWidget(centralWidget);

    // --- Connect Buttons ---
    connect(addBtn, &QPushButton::clicked, this, &MainWindow::addSequence);
    connect(searchBtn, &QPushButton::clicked, this, &MainWindow::searchKmer);
    connect(predictBtn, &QPushButton::clicked, this, &MainWindow::predictFunction);
    connect(gcBtn, &QPushButton::clicked, this, &MainWindow::showGCContent);
}

MainWindow::~MainWindow() {}

void MainWindow::addSequence() {
    QString fileName = QFileDialog::getOpenFileName(this, "Open FASTA File", "", "FASTA Files (*.fasta *.fa *.txt)");
    if (fileName.isEmpty())
        return;

    QFile file(fileName);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QMessageBox::warning(this, "Error", "Could not open file!");
        return;
    }

    QTextStream in(&file);
    QString content;
    while (!in.atEnd()) {
        QString line = in.readLine();
        if (!line.startsWith(">"))
            content += line.trimmed();
    }

    currentSequence = content.toUpper();
    resultWindow->setText("✅ Sequence loaded successfully!\n\nLength: " + QString::number(currentSequence.length()));
}

void MainWindow::searchKmer() {
    if (currentSequence.isEmpty()) {
        QMessageBox::warning(this, "Error", "No sequence loaded!");
        return;
    }
    bool ok;
    QString kmer = QInputDialog::getText(this, "Search K-mer", "Enter K-mer (any length)", QLineEdit::Normal, "", &ok);
    if (!ok || kmer.isEmpty())
        return;

    kmer = kmer.toUpper();
    int k = kmer.length();
    std::string genome = currentSequence.toStdString();
    std::string kmer_query = kmer.toStdString();
    KMerHashTable table(std::max(512, static_cast<int>(genome.length() / 2)), k);
    table.build_index(genome);

    std::vector<int> positions = table.query(kmer_query);
    QString output;
    if (positions.empty()) {
        output = QString("K-mer '%1' not found in the sequence.").arg(kmer);
    } else {
        output = QString("K-mer '%1' found %2 times:\n").arg(kmer).arg(positions.size());
        for (int p : positions)
            output += QString::number(p) + "\n";
    }
    // Optionally, show embedding prediction for this k-mer region
    resultWindow->setText(output);
}

double MainWindow::calculateGC(const QString &sequence) {
    int gcCount = 0;
    for (QChar c : sequence) {
        if (c == 'G' || c == 'C')
            gcCount++;
    }
    return (sequence.isEmpty()) ? 0.0 : (100.0 * gcCount / sequence.length());
}

void MainWindow::predictFunction() {
    if (currentSequence.isEmpty()) {
        QMessageBox::warning(this, "Error", "No sequence loaded!");
        return;
    }
    std::string genome = currentSequence.toStdString();
    int k = 5; // Example: use 5-mer embeddings
    KMerHashTable table(std::max(512, static_cast<int>(genome.length() / 2)), k);
    table.build_index(genome);
    auto counts = table.get_kmer_counts();
    SequenceEmbedding emb = get_sequence_embedding(counts, genome);

    double gc = calculateGC(currentSequence); // optional context
    std::string prediction = predict_function_from_embedding(emb, gc);

    resultWindow->setText(QString::fromStdString(prediction));
}

void MainWindow::showGCContent() {
    if (currentSequence.isEmpty()) {
        QMessageBox::warning(this, "Error", "No sequence loaded!");
        return;
    }

    double gc = calculateGC(currentSequence);
    resultWindow->setText(QString("📊 GC Content: %1%").arg(QString::number(gc, 'f', 2)));
}

