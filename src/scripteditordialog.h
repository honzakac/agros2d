#ifndef SCRIPTEDITORDIALOG_H
#define SCRIPTEDITORDIALOG_H

#include "scene.h"
#include "sceneview.h"
#include "scenemarker.h"
#include "scripteditorhighlighter.h"

class SceneView;

struct ScriptResult
{
    ScriptResult()
    {
        text = "";
        value = 0.0;
        isError = false;
    }

    ScriptResult(const QString &text, double value = 0, bool isError = false)
    {
        this->text = text;
        this->value = value;
        this->isError = isError;
    }

    QString text;
    double value;
    bool isError;
};

class Scene;
class SceneView;
class ScriptEditor;
class SearchDialog;

void createScripEngine(SceneView *sceneView);
QString createPythonFromModel();
ScriptResult runPython(const QString &script, const QString &fileName = "");
ScriptResult runPythonExpression(const QString &expression);

class ScriptEngineRemote : QObject
{
    Q_OBJECT

public:
    ScriptEngineRemote();
    ~ScriptEngineRemote();

private slots:
    void connected();
    void readCommand();
    void disconnected();

    void displayError(QLocalSocket::LocalSocketError socketError);

private:
    QString command;    

    QLocalServer *m_server;
    QLocalSocket *m_server_socket;
    QLocalSocket *m_client_socket;
};

class ScriptEditorWidget : public QWidget
{
    Q_OBJECT

public slots:
    void doRunPython();
    void doCreatePythonFromModel();
    void doPrintStdout(const QString &message);

public:
    QString file;
    ScriptEditor *txtEditor;
    QPlainTextEdit *txtOutput;

    ScriptEditorWidget(QWidget *parent);
    ~ScriptEditorWidget();

    void createControls();
    void createEngine();

private:
    QSplitter *splitter;
};

class ScriptEditorDialog : public QMainWindow
{
    Q_OBJECT
public:
    ScriptEditorDialog(QWidget *parent = 0);
    ~ScriptEditorDialog();

    void runScript(const QString &fileName);
    void runCommand(const QString &command);
    void showDialog();

public slots:
    void doFileNew();
    void doFileOpen(const QString &file = QString());
    void doFileSave();
    void doFileSaveAs();
    void doFileClose();
    void doFileOpenRecent(QAction *action);

    void doFind();
    void doFindNext(bool fromBegining = false);
    void doReplace();

    void doDataChanged();
    void doHelp();

    void doCloseTab(int index);
    void doCurrentPageChanged(int index);

private:
    QStringList recentFiles;

    ScriptEditor *txtEditor;
    SearchDialog *searchDialog;

    QMenu *mnuFile;
    QMenu *mnuRecentFiles;
    QMenu *mnuEdit;
    QMenu *mnuTools;
    QMenu *mnuHelp;

    QToolBar *tlbFile;
    QToolBar *tlbEdit;
    QToolBar *tlbTools;

    QAction *actFileNew;
    QAction *actFileOpen;
    QAction *actFileSave;
    QAction *actFileSaveAs;
    QAction *actFileClose;
    QAction *actExit;
    QActionGroup *actFileOpenRecentGroup;

    QAction *actUndo;
    QAction *actRedo;
    QAction *actCut;
    QAction *actCopy;
    QAction *actPaste;

    QAction *actFind;
    QAction *actFindNext;
    QAction *actReplace;

    QAction *actRunPython;
    QAction *actCreateFromModel;

    QAction *actHelp;

    QTabWidget *tabWidget;

    void createActions();
    void createControls();

    void setRecentFiles();
};

class ScriptEditor : public QPlainTextEdit
{
    Q_OBJECT

public:
    ScriptEditor(QWidget *parent = 0);
    ~ScriptEditor();

    void lineNumberAreaPaintEvent(QPaintEvent *event);
    int lineNumberAreaWidth();

protected:
    void resizeEvent(QResizeEvent *event);

private slots:
    void doUpdateLineNumberAreaWidth(int newBlockCount);
    void doHighlightCurrentLine();
    void doUpdateLineNumberArea(const QRect &, int);

private:
    QWidget *lineNumberArea;
};

class ScriptEditorLineNumberArea : public QWidget
{
public:
    ScriptEditorLineNumberArea(ScriptEditor *editor) : QWidget(editor) {
        codeEditor = editor;
    }

    QSize sizeHint() const {
        return QSize(codeEditor->lineNumberAreaWidth(), 0);
    }

protected:
    void paintEvent(QPaintEvent *event) {
        codeEditor->lineNumberAreaPaintEvent(event);
    }

private:
    ScriptEditor *codeEditor;
};

// ************************************************************************************************************

class SearchDialog: public QDialog
{
    Q_OBJECT
public:
    // Constructor
    SearchDialog(QWidget *parent=0);
    ~SearchDialog();

    int showDialogFind();
    int showDialogReplace();

    inline QString searchString()  { return txtFind->text(); }
    inline QString replaceString() { return txtReplace->text(); }

    inline bool searchStringIsRegExp() { return chkSearchRegExp->checkState(); }
    inline bool caseSensitive() { return chkCaseSensitive->checkState(); }

private:
    QLineEdit *txtFind, *txtReplace;
    QCheckBox *chkSearchRegExp, *chkCaseSensitive;
    QPushButton *btnCancel, *btnConfirm;
};

#endif // SCRIPTEDITORDIALOG_H
