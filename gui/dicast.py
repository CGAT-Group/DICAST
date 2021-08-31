import json
import os
from datetime import datetime
from pathlib import Path
from shutil import move, rmtree

import psutil
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import pyqtSignal, QTimer
from PyQt5.QtGui import QPixmap, QColor, QIntValidator
from PyQt5.QtWidgets import QAbstractItemView, QTabWidget, QWidget, QVBoxLayout, QPlainTextEdit

from dicast_config_setter import ConfigSetter, on_button_get_folder_path_clicked, on_button_get_file_path_clicked, show_popup

AS_TOOLS = ['asgal', 'aspli', 'eventpointer', 'irfinder', 'majiq', 'sgseq', 'spladder', 'whippet']
MAPPING_TOOLS = ['bbmap', 'contextmap', 'crac', 'dart', 'gsnap', 'hisat', 'mapsplice', 'minimap', 'segemehl', 'star', 'subjunc']


class Ui_main_window(object):
    trigger = pyqtSignal()

    def __init__(self, window):
        super().__init__()

        self.config_setter = None
        self.main_window = window
        self.main_grid = QtWidgets.QWidget(self.main_window)

        self.process = QtCore.QProcess()
        self._pid = -1
        self.process_finished = False
        self.std_err_file = ""
        self.std_out_file = ""
        self.new_process_starting_date = ""
        self.old_process_starting_date = ""

        self.timer = QTimer()

        self.resources_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), "resources")
        self.log_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), "log")

        self.statusbar = QtWidgets.QStatusBar(self.main_window)
        self.menubar = QtWidgets.QMenuBar(self.main_window)
        self.header_label = QtWidgets.QLabel(self.main_grid)

        self.line = QtWidgets.QFrame(self.main_grid)
        self.base_grid = QtWidgets.QGridLayout(self.main_grid)

        # File selection grid
        self.file_selector_grid = QtWidgets.QGridLayout()
        self.snake_label = QtWidgets.QLabel(self.main_grid)
        self.wd_label = QtWidgets.QLabel(self.main_grid)

        self.custom_snakefile_full_path = ""
        self.wd_folder_full_path = ""

        self.snake_button = QtWidgets.QPushButton(self.main_grid)
        self.snake_button.clicked.connect(lambda: self.set_snake_file(self.select_snakefile_label))

        self.wd_button = QtWidgets.QPushButton(self.main_grid)
        self.wd_button.clicked.connect(lambda: self.set_wd(self.select_wdfolder_label))

        self.select_snakefile_label = QtWidgets.QLabel(self.main_grid)
        self.select_wdfolder_label = QtWidgets.QLabel(self.main_grid)

        # Attribute grid
        self.attribute_settings_grid = QtWidgets.QGridLayout()
        self.set_configs_label = QtWidgets.QLabel(self.main_grid)
        self.set_configs_button = QtWidgets.QPushButton(self.main_grid)
        self.set_configs_button.clicked.connect(self.open_config_setter)
        self.config_setter_window = QtWidgets.QMainWindow()

        self.overwrite_label = QtWidgets.QLabel(self.main_grid)
        self.overwrite_cb = QtWidgets.QCheckBox(self.main_grid)
        self.overwrite_cb.stateChanged.connect(self.check_runnable)

        self.cores_label = QtWidgets.QLabel(self.main_grid)
        self.cores_line = QtWidgets.QLineEdit(self.main_grid)
        self.cores_line.setFixedWidth(120)
        self.cores_line.setText("2")
        self.cores_line.setValidator(QIntValidator())

        self.asim_label = QtWidgets.QLabel(self.main_grid)
        self.asim_cb = QtWidgets.QCheckBox(self.main_grid)

        self.mt_label = QtWidgets.QLabel(self.main_grid)
        self.mt_list = QtWidgets.QListWidget(self.main_grid)
        self.mt_list.setFixedHeight(130)
        self.mt_list.setSelectionMode(QAbstractItemView.MultiSelection)

        self.ast_label = QtWidgets.QLabel(self.main_grid)
        self.ast_list = QtWidgets.QListWidget(self.main_grid)
        self.ast_list.setFixedHeight(130)
        self.ast_list.setSelectionMode(QAbstractItemView.MultiSelection)

        for tool in MAPPING_TOOLS:
            self.mt_list.addItem(tool)

        for tool in AS_TOOLS:
            self.ast_list.addItem(tool)

        self.current_run_label = QtWidgets.QLabel(self.main_grid)
        self.refresh_status_button = QtWidgets.QPushButton(self.main_grid)
        self.refresh_status_button.clicked.connect(self.check_running)

        self.cleanup_button = QtWidgets.QPushButton(self.main_grid)
        self.cleanup_button.clicked.connect(self.cleanup)

        self.snake_out_text_edit = QPlainTextEdit()
        self.error_text_edit = QPlainTextEdit()
        self.snake_out_text_edit.setReadOnly(True)
        self.error_text_edit.setReadOnly(True)

        self.tabs = QTabWidget()
        self.tabs.setMinimumHeight(200)
        self.tabs.blockSignals(True)
        self.tabs.currentChanged.connect(self.refresh_log_display)
        self.snakemake_out_tab = QWidget()
        self.error_out_tab = QWidget()
        self.tabs.addTab(self.snakemake_out_tab, "DICAST Output")
        self.tabs.addTab(self.error_out_tab, "Snakemake Output")

        # Bottom button grid
        self.horizontal_layout = QtWidgets.QHBoxLayout()

        self.okcancel_buttonbox = QtWidgets.QDialogButtonBox(self.main_grid)
        self.okcancel_buttonbox.accepted.connect(self.on_button_ok_clicked)

    def setupUi(self, window):
        window.setObjectName("main_window")
        size_policy_f_f = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        size_policy_f_f.setHorizontalStretch(0)
        size_policy_f_f.setVerticalStretch(0)
        size_policy_e_f = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        size_policy_e_p = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        size_policy_e_p.setHorizontalStretch(0)
        size_policy_e_p.setVerticalStretch(0)

        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(window.sizePolicy().hasHeightForWidth())
        window.setSizePolicy(sizePolicy)
        icon = QtGui.QIcon()
        pix_map_path = os.path.join(self.resources_path, "Dicast.svg")
        icon.addPixmap(QtGui.QPixmap(pix_map_path), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        window.setWindowIcon(icon)
        self.main_grid.setEnabled(True)
        self.main_grid.setObjectName("main_grid")
        self.base_grid.setObjectName("gridLayout_2")
        horizontal_spacer = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.base_grid.addItem(horizontal_spacer, 6, 0, 1, 1)
        self.file_selector_grid.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.file_selector_grid.setObjectName("file_selector_grid")
        vertical_spacer = QtWidgets.QSpacerItem(0, 0, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)

        size_policy_e_p.setHeightForWidth(self.wd_label.sizePolicy().hasHeightForWidth())
        self.wd_label.setSizePolicy(size_policy_e_p)
        self.wd_label.setObjectName("wd_label")
        self.file_selector_grid.addWidget(self.wd_label, 0, 0, 1, 2, QtCore.Qt.AlignLeft)
        size_policy_f_f.setHeightForWidth(self.wd_button.sizePolicy().hasHeightForWidth())
        self.wd_button.setSizePolicy(size_policy_f_f)
        self.wd_button.setObjectName("wd_button")
        self.file_selector_grid.addWidget(self.wd_button, 1, 0, 1, 1, QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.file_selector_grid.addWidget(self.select_wdfolder_label, 1, 1, 1, 1, QtCore.Qt.AlignLeft | QtCore.Qt.AlignHCenter)
        self.file_selector_grid.addItem(vertical_spacer, 1, 3, 1, 1)

        size_policy_e_p.setHeightForWidth(self.snake_label.sizePolicy().hasHeightForWidth())
        self.snake_label.setSizePolicy(size_policy_e_p)
        self.snake_label.setObjectName("snake_label")
        self.file_selector_grid.addWidget(self.snake_label, 0, 2, 1, 2, QtCore.Qt.AlignLeft)
        size_policy_f_f.setHeightForWidth(self.snake_button.sizePolicy().hasHeightForWidth())
        self.snake_button.setSizePolicy(size_policy_f_f)
        self.snake_button.setObjectName("snake_button")
        self.snake_button.setToolTip("If you have a customized Snakefile, you can select it here."
                                     "\nOr leave it blank to use the one provided by DICAST."
                                     "\nIf you do use your own, the configuration file created with this"
                                     "\ntool will be saved to the same directory as your Snakefile.")
        self.file_selector_grid.addWidget(self.snake_button, 1, 2, 1, 1, QtCore.Qt.AlignLeft | QtCore.Qt.AlignHCenter)
        self.file_selector_grid.addWidget(self.select_snakefile_label, 1, 3, 1, 1, QtCore.Qt.AlignLeft | QtCore.Qt.AlignHCenter)
        self.file_selector_grid.addItem(vertical_spacer, 1, 1, 1, 1)

        size_policy_e_p.setHeightForWidth(self.set_configs_label.sizePolicy().hasHeightForWidth())
        self.set_configs_label.setSizePolicy(size_policy_e_p)
        self.set_configs_label.setObjectName("config_label")
        self.file_selector_grid.addWidget(self.set_configs_label, 0, 4, 1, 2, QtCore.Qt.AlignLeft)

        size_policy_f_f.setHeightForWidth(self.set_configs_button.sizePolicy().hasHeightForWidth())
        self.set_configs_button.setSizePolicy(size_policy_f_f)
        self.set_configs_button.setObjectName("config_button")
        self.file_selector_grid.addWidget(self.set_configs_button, 1, 4, 1, 1, QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        self.file_selector_grid.addItem(vertical_spacer, 1, 5, 1, 1)

        self.base_grid.addLayout(self.file_selector_grid, 2, 0, 1, 1)

        # Attributes
        # self.set_configs_label.setObjectName("set_configs_label")
        # self.attribute_settings_grid.addWidget(self.set_configs_label, 0, 0, 1, 1)

        # size_policy_f_f.setHeightForWidth(self.set_configs_button.sizePolicy().hasHeightForWidth())
        # self.set_configs_button.setSizePolicy(size_policy_f_f)
        # self.set_configs_button.setObjectName("set_configs_button")
        # self.attribute_settings_grid.addWidget(self.set_configs_button, 0, 1, 1, 1)

        self.overwrite_label.setObjectName("overwrite_label")
        self.attribute_settings_grid.addWidget(self.overwrite_label, 0, 0, 1, 1)

        self.overwrite_cb.setObjectName("overwrite_cb")
        self.attribute_settings_grid.addWidget(self.overwrite_cb, 0, 1, 1, 1, QtCore.Qt.AlignLeft)

        self.cores_label.setObjectName("cores_label")
        self.attribute_settings_grid.addWidget(self.cores_label, 1, 0, 1, 1)
        size_policy_e_f.setHeightForWidth(self.cores_line.sizePolicy().hasHeightForWidth())
        self.cores_line.setSizePolicy(size_policy_e_f)
        self.cores_line.setObjectName("cores_line")
        self.attribute_settings_grid.addWidget(self.cores_line, 1, 1, 1, 2)

        self.asim_label.setObjectName("asim_label")
        self.attribute_settings_grid.addWidget(self.asim_label, 2, 0, 1, 1)

        self.asim_cb.setObjectName("asim_cb")
        self.attribute_settings_grid.addWidget(self.asim_cb, 2, 1, 1, 1, QtCore.Qt.AlignLeft)

        self.attribute_settings_grid.setObjectName("attribute_settings_grid")
        self.mt_label.setObjectName("mt_label")
        self.attribute_settings_grid.addWidget(self.mt_label, 3, 0, 1, 1)

        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)

        self.mt_list.setMinimumHeight(60)
        self.mt_list.setSizePolicy(sizePolicy)
        self.mt_list.setObjectName("mt_list")
        self.attribute_settings_grid.addWidget(self.mt_list, 4, 0, 1, 1, QtCore.Qt.AlignLeft)

        self.ast_label.setObjectName("ast_label")
        self.attribute_settings_grid.addWidget(self.ast_label, 3, 1, 1, 1)

        self.ast_list.setMinimumHeight(60)
        self.ast_list.setSizePolicy(sizePolicy)
        self.ast_list.setObjectName("ast_list")
        self.attribute_settings_grid.addWidget(self.ast_list, 4, 1, 1, 1, QtCore.Qt.AlignLeft)

        self.current_run_label.setObjectName("current_run_label")
        self.current_run_label.setWordWrap(True)

        self.base_grid.addWidget(self.current_run_label, 5, 0, 1, 1)
        self.base_grid.addWidget(self.tabs, 6, 0, 1, 1)

        self.snake_out_text_edit.setObjectName("snake_out_text_edit")
        self.snakemake_out_tab.layout = QVBoxLayout()
        self.snakemake_out_tab.layout.addWidget(self.snake_out_text_edit, 0)
        self.snakemake_out_tab.setLayout(self.snakemake_out_tab.layout)

        self.error_text_edit.setObjectName("error_text_edit")
        self.error_out_tab.layout = QVBoxLayout()
        self.error_out_tab.layout.addWidget(self.error_text_edit, 0)
        self.error_out_tab.setLayout(self.error_out_tab.layout)

        self.horizontal_layout.setObjectName("horizontal_layout")

        size_policy_f_f.setHeightForWidth(self.refresh_status_button.sizePolicy().hasHeightForWidth())
        self.refresh_status_button.setSizePolicy(size_policy_f_f)
        self.refresh_status_button.setObjectName("refresh_status_button")
        self.horizontal_layout.addWidget(self.refresh_status_button)

        size_policy_f_f.setHeightForWidth(self.cleanup_button.sizePolicy().hasHeightForWidth())
        self.cleanup_button.setSizePolicy(size_policy_f_f)
        self.cleanup_button.setObjectName("refresh_status_button")
        self.cleanup_button.setToolTip("Should a run fail or be aborted, it may be necessary to remove\n"
                                       "unfinished log files and the lock on the working directory.")
        self.horizontal_layout.addWidget(self.cleanup_button)
        self.horizontal_layout.addItem(vertical_spacer)
        self.base_grid.addLayout(self.horizontal_layout, 7, 0, 1, 1)

        self.okcancel_buttonbox.setStandardButtons(QtWidgets.QDialogButtonBox.Close | QtWidgets.QDialogButtonBox.Abort | QtWidgets.QDialogButtonBox.Ok)
        self.submit_button = self.okcancel_buttonbox.button(QtWidgets.QDialogButtonBox.Ok)
        self.abort_button = self.okcancel_buttonbox.button(QtWidgets.QDialogButtonBox.Abort)
        palette = self.abort_button.palette()
        palette.setColor(self.abort_button.backgroundRole(), QColor(252, 53, 53))
        self.abort_button.setAutoFillBackground(True)
        self.abort_button.setPalette(palette)
        self.cancel_button = self.okcancel_buttonbox.button(QtWidgets.QDialogButtonBox.Close)
        self.abort_button.clicked.connect(self.on_button_abort_clicked)
        self.cancel_button.clicked.connect(self.main_window.close)
        self.okcancel_buttonbox.setObjectName("okcancel_buttonbox")
        self.base_grid.addWidget(self.okcancel_buttonbox, 8, 0, 1, 1)

        self.base_grid.addLayout(self.attribute_settings_grid, 3, 0, 1, 1)
        self.header_label.setObjectName("header_label")
        pix_map_path = os.path.join(self.resources_path, "Dicast-text.png")
        pixmap = QPixmap(pix_map_path)
        pixmap = pixmap.scaled(250, 250, QtCore.Qt.KeepAspectRatio)
        self.header_label.setPixmap(pixmap)
        self.base_grid.addWidget(self.header_label, 0, 0, 1, 1)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.base_grid.addWidget(self.line, 1, 0, 1, 1)
        window.setCentralWidget(self.main_grid)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 692, 20))
        self.menubar.setObjectName("menubar")
        window.setMenuBar(self.menubar)
        self.statusbar.setObjectName("statusbar")
        window.setStatusBar(self.statusbar)

        self.retranslateUi(window)
        QtCore.QMetaObject.connectSlotsByName(window)
        self.tabs.blockSignals(False)
        self.check_running()
        self.timer.timeout.connect(self.check_running)
        self.timer.setInterval(1000)

    def retranslateUi(self, window):
        _translate = QtCore.QCoreApplication.translate
        window.setWindowTitle(_translate("main_window", "DICAST"))
        self.wd_button.setText(_translate("main_window", "Browse"))
        self.snake_button.setText(_translate("main_window", "Browse"))
        self.snake_label.setText(_translate("main_window", "Select custom Snakefile"))
        self.select_snakefile_label.setText(_translate("main_window", "Select a file..."))
        self.select_wdfolder_label.setText(_translate("main_window", "Select a folder..."))
        self.wd_label.setText(_translate("main_window", "Select working directory"))
        self.mt_label.setText(_translate("main_window", "Which mapping tools: "))
        self.set_configs_label.setText(_translate("main_window", "Customize configuration file"))
        self.set_configs_button.setText(_translate("main_window", "Customize"))
        self.overwrite_label.setText(_translate("main_window", "Possible overwrite acknowledgment"))
        self.overwrite_cb.setText(_translate("main_window", "Acknowledged"))
        self.cores_label.setText(_translate("main_window", "Number of cores available to Snakemake"))
        self.ast_label.setText(_translate("main_window", "Which AS tools: "))
        self.asim_label.setText(_translate("main_window", "Do you want to run ASimulatoR?"))
        self.asim_cb.setText(_translate("main_window", "Yes"))
        self.refresh_status_button.setText(_translate("main_window", "Refresh status"))
        self.cleanup_button.setText(_translate("main_window", "Cleanup"))

    def set_wd(self, label):
        path = on_button_get_folder_path_clicked(label, "/")
        if path != "":
            self.wd_folder_full_path = path
        wd_button_palette = self.wd_button.style().standardPalette()
        if self.wd_folder_full_path != "":
            self.config_setter = ConfigSetter(self.config_setter_window, self.wd_folder_full_path)
            if not self.config_setter.parse_successful:
                self.wd_folder_full_path = ""
                self.select_wdfolder_label.setText("Select a folder...")
                self.config_setter = None
                wd_button_palette.setColor(self.wd_button.backgroundRole(), QColor(252, 53, 53))
            else:
                self.config_setter.setup_ui(self.config_setter_window)
            self.wd_button.setPalette(wd_button_palette)

    def set_snake_file(self, label):
        self.custom_snakefile_full_path = on_button_get_file_path_clicked(label)

    def on_button_ok_clicked(self):
        if self.check_runnable(True):
            self.submit()
            self.check_running()

    def submit(self):
        self.start_process()
        self.process.started.connect(lambda: self.submit_button.setDisabled(True))
        self.process.finished.connect(lambda: self.submit_button.setEnabled(True))
        self.process.started.connect(lambda: self.cleanup_button.setDisabled(True))
        self.process.finished.connect(lambda: self.cleanup_button.setEnabled(True))

    def check_runnable(self, submit_state=False):
        overwrite_cb_palette = self.overwrite_cb.palette()
        if submit_state:
            if not self.overwrite_cb.isChecked():
                overwrite_cb_palette.setColor(self.overwrite_cb.backgroundRole(), QColor(252, 128, 121))
                self.overwrite_cb.setAutoFillBackground(True)
            else:
                overwrite_cb_palette = self.overwrite_cb.style().standardPalette()

            wd_button_palette = self.wd_button.palette()
            if self.wd_folder_full_path == "":
                wd_button_palette.setColor(self.wd_button.backgroundRole(), QColor(252, 53, 53))
            else:
                wd_button_palette = self.wd_button.style().standardPalette()
            self.wd_button.setPalette(wd_button_palette)

        self.overwrite_cb.setPalette(overwrite_cb_palette)

        return self.overwrite_cb.isChecked()

    def check_running(self):
        def check_process():
            log_file = os.path.join(self.log_path, "_pid.txt")
            try:
                with open(log_file, "r") as f:
                    process_identification = f.read()
                    _pid = int(process_identification.split("-")[0])
                    self.old_process_starting_date = process_identification.split("-")[1]
                    try:
                        if _pid != "" and _pid != -1:
                            self._pid = _pid
                            process = psutil.Process(_pid)

                            # Check if the process with the PID parsed from the _pid.txt file is actually the DICAST process by comparing it's creation time
                            # with the datetime parsed from the _pid.txt file. If the difference is smaller than 5 seconds, we assume its the same process
                            # 5 seconds because there can be a time delay between setting the datetime for the log files and the creation of the process
                            old_datetime = datetime.strptime(self.old_process_starting_date, '%Y_%m_%d_%H_%M_%S')
                            time_difference = abs((old_datetime - datetime.fromtimestamp(process.create_time())).total_seconds())
                            if time_difference < 5:
                                p_status = psutil.Process(_pid).status()
                                self.submit_button.setEnabled(False)
                            else:
                                raise psutil.NoSuchProcess
                        else:
                            return "Current status: Empty _pid.txt"
                    except psutil.NoSuchProcess:
                        self.submit_button.setEnabled(True)
                        if self.timer.isActive():
                            self.timer.stop()
                        return "Current status: No instance running, previous instance finished. " \
                               "Tabs below will show log files from the previous run."
                    else:
                        if not self.timer.isActive():
                            self.timer.start()
                        return f"Current status: Instance {_pid} has status: {p_status}"
            except FileNotFoundError:
                return f"Current status: Missing pid file. No instance running or wrong folder: {log_file}"
        self.current_run_label.setText(check_process())
        self.run_all_refresh_log_display()

    def run_all_refresh_log_display(self):
        self.refresh_log_display(0)
        self.refresh_log_display(1)

    def build_snakemake_config(self):
        config_dict = {
            "Possible_overwrite_acknowledge": self.overwrite_cb.isChecked(),
            "ASimulatoR": self.asim_cb.isChecked(),
            "Mapping_tools": " ".join([item.text() for item in self.mt_list.selectedItems()]),
            "Alternative_splicing_detection_tools": " ".join([item.text() for item in self.ast_list.selectedItems()])
        }

        if self.custom_snakefile_full_path != "":
            config_json_path = os.path.dirname(self.custom_snakefile_full_path)
            config_json_path = os.path.join(config_json_path, 'snakemake_config.json')
        else:
            config_json_path = os.path.join(self.wd_folder_full_path, 'scripts', 'snakemake', 'snakemake_config.json')
            if Path(config_json_path).is_file():
                date = datetime.today().strftime('%Y_%m_%d_%H_%M_%S')
                backup_config_file = os.path.join(os.path.dirname(config_json_path), f"backup_{date}_{os.path.basename(config_json_path)}")
                os.rename(config_json_path, backup_config_file)
                backup_path = os.path.join(os.path.dirname(config_json_path), 'backup')
                if not os.path.exists(backup_path):
                    os.makedirs(backup_path)
                move(backup_config_file, backup_path)

        with open(config_json_path, 'w') as f:
            json.dump(config_dict, f)

    def create_command_line_args(self):
        if self.custom_snakefile_full_path != "":
            snakefile_path = self.custom_snakefile_full_path
            snakefile_config_path = os.path.dirname(self.custom_snakefile_full_path)
        else:
            snakefile_path = os.path.join(self.wd_folder_full_path, 'scripts', 'snakemake', 'Snakefile')
            snakefile_config_path = os.path.join(self.wd_folder_full_path, 'scripts', 'snakemake', 'snakemake_config.json')

        if not Path(snakefile_path).is_file():
            show_popup("Error - File not found", f"Couldn't find Snakefile {snakefile_path}."
                                                 f"\nPlease make sure you have selected the correct working directory and the file is in the correct folder.")
            return False

        if not Path(snakefile_config_path).is_file():
            show_popup("Error - File not found", f"Couldn't find Snakefile config {snakefile_config_path}."
                                                 f"\nSomething must have gone wrong in the creation of the config file. Please make sure the working directory is correct."
                                                 f"\nIf you used a custom Snakefile, make sure you have write access to the folder.")
            return False

        args = ["-j", f"{self.cores_line.text()}", "-s", f"{snakefile_path}", "-d", f"{self.wd_folder_full_path}", "--configfile", f"{snakefile_config_path}"]
        return args

    def start_process(self):
        self.build_snakemake_config()
        run_str = "snakemake"
        args = self.create_command_line_args()
        if not args:
            return
        self.process.setProgram(run_str)
        self.process.setArguments(args)
        if not os.path.isdir(self.log_path):
            os.makedirs(self.log_path)
        self.new_process_starting_date = datetime.today().strftime('%Y_%m_%d_%H_%M_%S')
        self.std_err_file = os.path.join(self.log_path, f"snakemake_log_{self.new_process_starting_date}.txt")
        self.std_out_file = os.path.join(self.log_path, f"dicast_log_{self.new_process_starting_date}.txt")
        self.process.setStandardOutputFile(self.std_out_file)
        self.process.setStandardErrorFile(self.std_err_file)
        ok, pid = self.process.startDetached()
        if ok:
            self._pid = pid
            with open(os.path.join(self.log_path, "_pid.txt"), "w") as f:
                process_identification = f"{str(self._pid)}-{self.new_process_starting_date}"
                f.write(process_identification)
        else:
            show_popup("Error - Failed to start.", "An error occurred and DICAST was unable to run the Snakemake pipeline.\n"
                                                   "Make sure you have activated the correct conda environment and you can run snakemake.\n"
                                                   "See the documentation for more information.")

    def on_button_abort_clicked(self):
        if not self.process_finished:
            self.stop_process()
            self.check_running()

    def process_finished(self):
        self.submit_button.setEnabled(True)
        self.process_finished = True

    def stop_process(self):
        if self._pid > 0:
            try:
                p = psutil.Process(self._pid)
            except psutil.NoSuchProcess:
                pass
            else:
                for child in p.children(recursive=True):
                    child.kill()
                p.kill()
                self._pid = -1

    def cleanup(self):
        if self.wd_folder_full_path == "":
            return
        locks_folder = os.path.join(self.wd_folder_full_path, ".snakemake", "locks")
        snake_output_folder = os.path.join(self.wd_folder_full_path, "output", "snakemake")
        if os.path.exists(locks_folder):
            rmtree(locks_folder)
        if os.path.exists(snake_output_folder):
            rmtree(snake_output_folder)

    def refresh_log_display(self, i):
        if i == 0:
            try:
                with open(self.std_out_file, 'r') as dicast_log_file:
                    data = dicast_log_file.read()
                    self.snake_out_text_edit.setPlainText(data)
                    # Scroll to the bottom
                    self.snake_out_text_edit.verticalScrollBar().setValue(self.snake_out_text_edit.verticalScrollBar().maximum())
            except FileNotFoundError:
                try:
                    log_file = os.path.join(self.log_path, "dicast_log_" + self.old_process_starting_date + ".txt")
                    with open(log_file, 'r') as out_dicast_log_file:
                        data = out_dicast_log_file.read()
                        self.snake_out_text_edit.setPlainText(data)
                        self.snake_out_text_edit.verticalScrollBar().setValue(self.snake_out_text_edit.verticalScrollBar().maximum())
                except FileNotFoundError:
                    self.snake_out_text_edit.setPlainText("Couldn't find DICAST log file from previous run.")
        elif i == 1:
            try:
                with open(self.std_err_file, 'r') as out_snakemake_log_file:
                    data = out_snakemake_log_file.read()
                    self.error_text_edit.setPlainText(data)
                    self.error_text_edit.verticalScrollBar().setValue(self.error_text_edit.verticalScrollBar().maximum())

            except FileNotFoundError:
                try:
                    snakemake_log_file = os.path.join(self.log_path, "snakemake_log_" + self.old_process_starting_date + ".txt")
                    with open(snakemake_log_file, 'r') as out_snakemake_log_file:
                        data = out_snakemake_log_file.read()
                        self.error_text_edit.setPlainText(data)
                        self.error_text_edit.verticalScrollBar().setValue(self.error_text_edit.verticalScrollBar().maximum())
                except FileNotFoundError:
                    self.error_text_edit.setPlainText("Couldn't find Snakemake log file from previous run.")

    def open_config_setter(self):
        if self.config_setter is not None and self.config_setter.parse_successful:
            self.config_setter_window.show()
        else:
            show_popup("Error - No config file found", "No config file found. Make sure the correct working directory is set.")


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    main_window = QtWidgets.QMainWindow()
    main_window.resize(600, 800)
    ui = Ui_main_window(main_window)
    ui.setupUi(main_window)
    main_window.show()
    ui.run_all_refresh_log_display()
    sys.exit(app.exec_())
