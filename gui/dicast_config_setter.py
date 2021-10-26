import os
import re
from datetime import datetime
from pathlib import Path
from shutil import copymode, move
from tempfile import mkstemp
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtGui import QIntValidator, QPixmap, QColor
from PyQt5.QtWidgets import QMessageBox, QFileDialog, QToolButton, QPushButton


class FileDialog(QFileDialog):
    def __init__(self, restrictive):
        super().__init__()
        if restrictive:
            self.forward_button = self.findChild(QToolButton, "forwardButton")
            self.forward_button.blockSignals(True)
            self.back_button = self.findChild(QToolButton, "backButton")
            self.back_button.blockSignals(True)


def on_button_get_folder_path_clicked(label, path_restriction):
    if path_restriction == "/":
        dialog = FileDialog(False)
    else:
        dialog = FileDialog(True)
    dialog.setOption(QFileDialog.DontUseNativeDialog)
    dialog.setFileMode(QFileDialog.Directory)
    if path_restriction == "/":
        dialog.setDirectory(".")
    else:
        if os.path.exists(path_restriction):
            dialog.setDirectory(path_restriction)
        else:
            show_popup("Error - Folder not found", f'This folder has to be located in "{path_restriction}", however, this folder couldn\'t be found.')
            return ""
    dialog.directoryEntered.connect(lambda p: dir_changed(p, path_restriction, dialog))
    dialog.exec_()
    if len(dialog.selectedFiles()) == 0:
        return ""
    folder = dialog.selectedFiles()[0]
    if folder != "":
        label.setText(os.path.basename(folder))
        return folder
    return ""


def dir_changed(path, restriction, dialog):
    # https://stackoverflow.com/a/37095733
    def path_is_allowed(allowed_path, entered_path):
        allowed_path = os.path.abspath(allowed_path)
        entered_path = os.path.abspath(entered_path)
        common_path = os.path.commonpath([allowed_path, entered_path])
        if entered_path == allowed_path or common_path == allowed_path:
            return True
        else:
            return False

    if not path_is_allowed(restriction, path):
        dialog.setDirectory(restriction)


def on_button_get_file_path_clicked(label, path_restriction):
    if path_restriction == "/":
        dialog = FileDialog(False)
    else:
        dialog = FileDialog(True)
    dialog.setOption(QFileDialog.DontUseNativeDialog)
    dialog.setFileMode(QFileDialog.ExistingFile)
    if path_restriction == "/":
        dialog.setDirectory(".")
    else:
        if os.path.exists(path_restriction):
            dialog.setDirectory(path_restriction)
        else:
            show_popup("Error - Folder not found", f'This file has to be located in "{path_restriction}", however, this folder couldn\'t be found.')
            return ""
    dialog.directoryEntered.connect(lambda p: dir_changed(p, path_restriction, dialog))
    dialog.exec_()
    if len(dialog.selectedFiles()) == 0:
        return ""
    path = dialog.selectedFiles()[0]
    if path != "":
        label.setText(os.path.basename(path))
        return path
    return ""


def show_popup(title, message, error=True):
    msg = QMessageBox()
    msg.setWindowTitle(title)
    msg.setText(message)
    if error:
        msg.setIcon(QMessageBox.Critical)
    msg.exec_()


def on_button_close_clicked(window):
    mb = QMessageBox(QMessageBox.Question, "Confirm close", "Are you sure you want to close? All unsaved changes will be lost.")
    mb.addButton(QMessageBox.Yes)
    mb.addButton(QMessageBox.No)
    mb.setDefaultButton(QMessageBox.No)
    reply = mb.exec()
    if reply == QMessageBox.Yes:
        window.close()


class ConfigSetter(object):

    def __init__(self, window, working_directory):
        self.wd = working_directory
        self.parse_successful = False
        if not os.path.exists(os.path.join(self.wd, "input")):
            show_popup("Error - Missing input directory", f"Couldn't find directory {os.path.join(self.wd, 'input')}."
                                                          f"\nPlease make sure you have selected the correct working directory and the input directory is called 'input'.")
            return
        self.main_config_file = ""
        self.as_config_file = ""
        self.main_config_dict = {}
        self.as_config_dict = {}
        self.main_window = window

        self.resources_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), "resources")

        self.parse_successful = self.can_open_config_files()
        if not self.parse_successful:
            return
        self.centralwidget = QtWidgets.QWidget(self.main_window)
        self.header_label = QtWidgets.QLabel(self.centralwidget)
        self.line = QtWidgets.QFrame(self.centralwidget)

        self.main_window.resize(450, 450)
        self.centralwidget.setEnabled(True)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.centralwidget)
        self.gridLayout = QtWidgets.QGridLayout()
        self.onlyInt = QIntValidator()

        self.ncores = QtWidgets.QLabel(self.centralwidget)
        self.ncores_line = QtWidgets.QLineEdit(self.centralwidget)
        self.ncores_line.setValidator(self.onlyInt)

        self.read_length = QtWidgets.QLabel(self.centralwidget)
        self.read_length_line = QtWidgets.QLineEdit(self.centralwidget)
        self.read_length_line.setValidator(self.onlyInt)

        self.out_dir = QtWidgets.QLabel(self.centralwidget)
        self.out_dir_button = QtWidgets.QPushButton(self.centralwidget)
        self.out_dir_path = QtWidgets.QLabel(self.centralwidget)
        self.out_dir_full_path = ""
        self.out_dir_button.clicked.connect(lambda: self.set_out_dir(self.out_dir_path))

        self.fasta_name = QtWidgets.QLabel(self.centralwidget)
        self.fasta_name_button = QtWidgets.QPushButton(self.centralwidget)
        self.fasta_name_path = QtWidgets.QLabel(self.centralwidget)
        self.fasta_full_path = ""
        self.fasta_name_button.clicked.connect(lambda: self.set_fasta(self.fasta_name_path))

        self.gtf_name = QtWidgets.QLabel(self.centralwidget)
        self.gtf_name_button = QtWidgets.QPushButton(self.centralwidget)
        self.gtf_name_path = QtWidgets.QLabel(self.centralwidget)
        self.gtf_full_path = ""
        self.gtf_name_button.clicked.connect(lambda: self.set_gtf(self.gtf_name_path))

        self.gff_name = QtWidgets.QLabel(self.centralwidget)
        self.gff_name_button = QtWidgets.QPushButton(self.centralwidget)
        self.gff_name_path = QtWidgets.QLabel(self.centralwidget)
        self.gff_full_path = ""
        self.gff_name_button.clicked.connect(lambda: self.set_gff(self.gff_name_path))

        self.fastq_pair1_suffix = QtWidgets.QLabel(self.centralwidget)
        self.fastq_pair1_suffix_line = QtWidgets.QLineEdit(self.centralwidget)
        self.fastq_pair2_suffix = QtWidgets.QLabel(self.centralwidget)
        self.fastq_pair2_suffix_line = QtWidgets.QLineEdit(self.centralwidget)

        self.use_bam_input_files = QtWidgets.QLabel(self.centralwidget)
        self.use_bam_input_files_cb = QtWidgets.QCheckBox(self.centralwidget)

        self.combine_events_cb = QtWidgets.QCheckBox(self.centralwidget)
        self.combine_events = QtWidgets.QLabel(self.centralwidget)

        self.save_cancel_button_box = QtWidgets.QDialogButtonBox(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(window)
        self.config_translate_dict = {}

    def setup_ui(self, window):
        window.setObjectName("config_setter")
        self.verticalLayout.setObjectName("verticalLayout")
        self.header_label.setObjectName("header_label")
        pix_map_path = os.path.join(self.resources_path, "dicast_setter_text.png")
        pixmap = QPixmap(pix_map_path)
        pixmap = pixmap.scaled(250, 250, QtCore.Qt.KeepAspectRatio)
        self.header_label.setPixmap(pixmap)
        self.verticalLayout.addWidget(self.header_label)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.verticalLayout.addWidget(self.line)

        icon = QtGui.QIcon()
        pix_map_path = os.path.join(self.resources_path, "dicast_setter.svg")
        icon.addPixmap(QtGui.QPixmap(pix_map_path), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        window.setWindowIcon(icon)

        size_policy_e_f = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        size_policy_max_f = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Fixed)
        size_policy_min_f = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)

        self.gridLayout.setHorizontalSpacing(10)
        self.gridLayout.setVerticalSpacing(12)
        self.gridLayout.setObjectName("gridLayout")

        self.ncores.setObjectName("ncores")
        self.gridLayout.addWidget(self.ncores, 0, 0, 1, 1)
        size_policy_e_f.setHeightForWidth(self.ncores_line.sizePolicy().hasHeightForWidth())
        self.ncores_line.setSizePolicy(size_policy_e_f)
        self.ncores_line.setObjectName("ncores_line")
        self.gridLayout.addWidget(self.ncores_line, 0, 1, 1, 2)

        self.read_length.setObjectName("read_length")
        self.gridLayout.addWidget(self.read_length, 1, 0, 1, 1)
        size_policy_e_f.setHeightForWidth(self.read_length_line.sizePolicy().hasHeightForWidth())
        self.read_length_line.setSizePolicy(size_policy_e_f)
        self.read_length_line.setObjectName("read_length_line")
        self.gridLayout.addWidget(self.read_length_line, 1, 1, 1, 2)

        self.out_dir.setObjectName("out_dir")
        self.gridLayout.addWidget(self.out_dir, 2, 0, 1, 1)
        size_policy_max_f.setHeightForWidth(self.out_dir_button.sizePolicy().hasHeightForWidth())
        self.out_dir_button.setSizePolicy(size_policy_max_f)
        self.out_dir_button.setObjectName("out_dir_button")
        self.gridLayout.addWidget(self.out_dir_button, 2, 1, 1, 1)
        self.out_dir_path.setObjectName("out_dir_path")
        self.gridLayout.addWidget(self.out_dir_path, 2, 2, 1, 1)

        self.fasta_name.setObjectName("fasta_name")
        self.gridLayout.addWidget(self.fasta_name, 3, 0, 1, 1)
        self.fasta_name_path.setObjectName("fasta_name_path")
        self.gridLayout.addWidget(self.fasta_name_path, 3, 2, 1, 1)
        size_policy_max_f.setHeightForWidth(self.fasta_name_button.sizePolicy().hasHeightForWidth())
        self.fasta_name_button.setSizePolicy(size_policy_max_f)
        self.fasta_name_button.setObjectName("fasta_name_button")
        self.gridLayout.addWidget(self.fasta_name_button, 3, 1, 1, 1)

        self.gtf_name.setObjectName("gtf_name")
        self.gridLayout.addWidget(self.gtf_name, 4, 0, 1, 1)
        self.gtf_name_path.setObjectName("gtf_name_path")
        self.gridLayout.addWidget(self.gtf_name_path, 4, 2, 1, 1)
        size_policy_max_f.setHeightForWidth(self.gtf_name_button.sizePolicy().hasHeightForWidth())
        self.gtf_name_button.setSizePolicy(size_policy_max_f)
        self.gtf_name_button.setObjectName("gtf_name_button")
        self.gridLayout.addWidget(self.gtf_name_button, 4, 1, 1, 1)

        self.gff_name.setObjectName("gff_name")
        self.gridLayout.addWidget(self.gff_name, 5, 0, 1, 1)
        self.gff_name_path.setObjectName("gff_name_path")
        self.gridLayout.addWidget(self.gff_name_path, 5, 2, 1, 1)
        size_policy_max_f.setHeightForWidth(self.gff_name_button.sizePolicy().hasHeightForWidth())
        self.gff_name_button.setSizePolicy(size_policy_max_f)
        self.gff_name_button.setObjectName("gff_name_button")
        self.gridLayout.addWidget(self.gff_name_button, 5, 1, 1, 1)

        self.fastq_pair1_suffix.setObjectName("fastq_pair1_suffix")
        self.gridLayout.addWidget(self.fastq_pair1_suffix, 6, 0, 1, 1)
        size_policy_e_f.setHeightForWidth(self.fastq_pair1_suffix_line.sizePolicy().hasHeightForWidth())
        self.fastq_pair1_suffix_line.setSizePolicy(size_policy_e_f)
        self.fastq_pair1_suffix_line.setObjectName("fastq_pair1_suffix_line")
        self.gridLayout.addWidget(self.fastq_pair1_suffix_line, 6, 1, 1, 2)

        self.fastq_pair2_suffix.setObjectName("fastq_pair2_suffix")
        self.gridLayout.addWidget(self.fastq_pair2_suffix, 7, 0, 1, 1)
        size_policy_e_f.setHeightForWidth(self.fastq_pair2_suffix_line.sizePolicy().hasHeightForWidth())
        self.fastq_pair2_suffix_line.setSizePolicy(size_policy_e_f)
        self.fastq_pair2_suffix_line.setObjectName("fastq_pair2_suffix_line")
        self.gridLayout.addWidget(self.fastq_pair2_suffix_line, 7, 1, 1, 2)

        self.use_bam_input_files.setObjectName("use_bam_input_files")
        self.gridLayout.addWidget(self.use_bam_input_files, 8, 0, 1, 1)
        size_policy_min_f.setHeightForWidth(self.use_bam_input_files_cb.sizePolicy().hasHeightForWidth())
        self.use_bam_input_files_cb.setSizePolicy(size_policy_min_f)
        self.use_bam_input_files_cb.setObjectName("use_bam_input_files_cb")
        self.gridLayout.addWidget(self.use_bam_input_files_cb, 8, 1, 1, 1)

        self.combine_events.setObjectName("combine_events")
        self.combine_events.setWordWrap(True)
        self.gridLayout.addWidget(self.combine_events, 9, 0, 1, 1)
        size_policy_min_f.setHeightForWidth(self.combine_events_cb.sizePolicy().hasHeightForWidth())
        self.combine_events_cb.setSizePolicy(size_policy_min_f)
        self.combine_events_cb.setObjectName("combine_events_cb")

        self.gridLayout.addWidget(self.combine_events_cb, 9, 1, 1, 1)
        vertical_spacer = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        horizontal_spacer = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)

        self.gridLayout.addItem(horizontal_spacer)
        self.verticalLayout.addLayout(self.gridLayout)
        self.verticalLayout.addItem(vertical_spacer)

        self.save_cancel_button_box.setStandardButtons(QtWidgets.QDialogButtonBox.Close | QtWidgets.QDialogButtonBox.Save)
        self.save_cancel_button_box.setObjectName("save_cancel_button_box")
        self.save_button = self.save_cancel_button_box.button(QtWidgets.QDialogButtonBox.Save)
        self.save_button.clicked.connect(self.on_button_save_clicked)
        self.close_button = self.save_cancel_button_box.button(QtWidgets.QDialogButtonBox.Close)
        self.close_button.clicked.connect(lambda: on_button_close_clicked(self.main_window))
        self.verticalLayout.addWidget(self.save_cancel_button_box)

        window.setCentralWidget(self.centralwidget)
        self.statusbar.setEnabled(True)
        self.statusbar.setObjectName("statusbar")
        window.setStatusBar(self.statusbar)

        self.parse_old_configs_to_dicts()
        self.retranslateUi(window)
        QtCore.QMetaObject.connectSlotsByName(window)

    def retranslateUi(self, window):
        _translate = QtCore.QCoreApplication.translate
        window.setWindowTitle(_translate("config_setter", "Configuration Setter"))
        self.fasta_name.setText(_translate("config_setter", "Fasta name"))
        self.fastq_pair1_suffix.setText(_translate("config_setter", "First fastq pair suffix"))
        self.fastq_pair2_suffix.setText(_translate("config_setter", "Second fastq pair suffix"))
        self.read_length.setText(_translate("config_setter", "Read length"))
        self.gtf_name_button.setText(_translate("config_setter", "Browse"))
        self.fasta_name_button.setText(_translate("config_setter", "Browse"))
        self.gtf_name.setText(_translate("config_setter", "GTF name"))
        self.use_bam_input_files_cb.setText(_translate("config_setter", "*.bam files"))
        self.use_bam_input_files.setText(_translate("config_setter", "Use BAM files as input"))
        self.combine_events_cb.setText(_translate("config_setter", "Combine"))
        self.out_dir_button.setText(_translate("config_setter", "Browse"))
        self.out_dir.setText(_translate("config_setter", "Base output directory"))
        self.gff_name_button.setText(_translate("config_setter", "Browse"))
        self.combine_events.setText(_translate("config_setter", "Should different AS events be combined into one"))
        self.ncores.setText(_translate("config_setter", "Number of cores for each tool"))
        self.gff_name.setText(_translate("config_setter", "GFF name"))

        # Show default values in text fields
        self.ncores_line.setText(self.main_config_dict["ncores"])
        self.read_length_line.setText(self.main_config_dict["read_length"])

        self.out_dir_full_path = self.main_config_dict["outdir"]
        self.fasta_full_path = self.main_config_dict["fastaname"]
        self.gtf_full_path = self.main_config_dict["gtfname"]
        self.gff_full_path = self.main_config_dict["gffname"]

        self.out_dir_path.setText(_translate("config_setter", os.path.basename(self.out_dir_full_path)))
        self.fasta_name_path.setText(_translate("config_setter", os.path.basename(self.fasta_full_path)))
        self.gtf_name_path.setText(_translate("config_setter", os.path.basename(self.gtf_full_path)))
        self.gff_name_path.setText(_translate("config_setter", os.path.basename(self.gff_full_path)))

        self.fastq_pair1_suffix_line.setText(self.as_config_dict["fastqpair1suffix"])
        self.fastq_pair2_suffix_line.setText(self.as_config_dict["fastqpair2suffix"])
        self.use_bam_input_files_cb.setChecked(True if self.as_config_dict["use_bam_input_files"] == "1" else False)
        self.combine_events_cb.setChecked(True if self.as_config_dict["combine_events"] == "1" else False)

    def can_open_config_files(self):
        self.main_config_file = os.path.join(self.wd, "scripts", "config.sh")
        self.as_config_file = os.path.join(self.wd, "scripts", "asevent_config.sh")

        if Path(self.main_config_file).is_file() and Path(self.as_config_file).is_file():
            return True

        show_popup("Error - File not found", f"Couldn't find config files {self.main_config_file} and {self.as_config_file}. "
                                             f"\nPlease make sure you have selected the correct working directory.")
        return False

    def update_translate_dict(self):
        self.config_translate_dict = {
            "ncores": self.ncores_line.text(),
            "read_length": self.read_length_line.text(),
            "outdir": self.out_dir_full_path,
            "fastaname": os.path.basename(self.fasta_full_path),
            "gtfname": os.path.basename(self.gtf_full_path),
            "gffname": os.path.basename(self.gff_full_path),
            "fastqpair1suffix": self.fastq_pair1_suffix_line.text(),
            "fastqpair2suffix": self.fastq_pair2_suffix_line.text(),
            "use_bam_input_files": "1" if self.use_bam_input_files_cb.isChecked() else "0",
            "combine_events": "1" if self.combine_events_cb.isChecked() else "0"
        }

    def on_button_save_clicked(self):
        try:
            self.update_translate_dict()
            with open(self.main_config_file, 'r') as original_main_config:
                with open(self.as_config_file, 'r') as original_as_config:
                    self.write_new_config_file(original_as_config, self.as_config_file)
                    self.write_new_config_file(original_main_config, self.main_config_file)
            return True
        except FileNotFoundError:
            show_popup("Error - File not found", f"Couldn't find config files {self.main_config_file} and {self.as_config_file}."
                                                 f"\nPlease make sure you have selected the correct working directory.")
            return False

    def parse_old_configs_to_dicts(self):
        def parse_old_config_to_dict(config, config_dict):
            for line in config:
                if '=' in line and not line.startswith('#'):
                    line_arr = line.split('=')
                    config_key = line_arr[0].strip()
                    regex = re.compile("^([^\s|\t]+)")
                    config_value = regex.match(line_arr[1]).group(1)
                    config_dict[config_key] = config_value

        try:
            with open(self.main_config_file, 'r') as original_main_config:
                with open(self.as_config_file, 'r') as original_as_config:
                    parse_old_config_to_dict(original_main_config, self.main_config_dict)
                    parse_old_config_to_dict(original_as_config, self.as_config_dict)
        except FileNotFoundError:
            show_popup("Error - File not found", f"Couldn't find config files {self.main_config_file} and {self.as_config_file}."
                                                 f"\nPlease make sure you have selected the correct working directory.")

    def write_new_config_file(self, original_config, original_config_path):
        new_config_file_handler, new_config_path = mkstemp()
        with os.fdopen(new_config_file_handler, 'w') as new_config:
            config_dict = {}
            for line in original_config:
                if '=' in line and not line.startswith('#'):
                    line_arr = line.split('=')
                    try:
                        config_key = line_arr[0].strip()
                        new_config_value = self.config_translate_dict[config_key]
                    except KeyError:
                        new_config.write(line)
                        continue
                    else:
                        if new_config_value == "":
                            new_config.write(line)
                            continue

                    regex = re.compile("^([^\s|\t]+)")
                    config_value = regex.match(line_arr[1]).group(1)
                    config_dict[config_key] = config_value
                    new_config_value_line = line.replace(config_value, new_config_value)
                    new_config.write(new_config_value_line)
                else:
                    new_config.write(line)
        copymode(original_config_path, new_config_path)
        date = datetime.today().strftime('%Y_%m_%d_%H_%M_%S')
        backup_config_file = os.path.join(os.path.dirname(original_config_path), f"backup_{date}_{os.path.basename(original_config_path)}")
        os.rename(original_config_path, backup_config_file)
        move(new_config_path, original_config_path)
        backup_path = os.path.join(os.path.dirname(backup_config_file), 'backup')
        if not os.path.exists(backup_path):
            os.makedirs(backup_path)
        move(backup_config_file, backup_path)

    def set_fasta(self, path):
        self.fasta_full_path = on_button_get_file_path_clicked(path, os.path.join(self.wd, "input"))

    def set_gtf(self, path):
        self.gtf_full_path = on_button_get_file_path_clicked(path, os.path.join(self.wd, "input"))

    def set_gff(self, path):
        self.gff_full_path = on_button_get_file_path_clicked(path, os.path.join(self.wd, "input"))

    def set_out_dir(self, path):
        self.out_dir_full_path = on_button_get_folder_path_clicked(path, self.wd)
        mount_path = self.out_dir_full_path.replace(self.wd, "$workdir")
        self.out_dir_full_path = os.path.join(mount_path, "${tool:-unspecific}-output")


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    config_setter = QtWidgets.QMainWindow()
    ui = ConfigSetter(config_setter, str(sys.argv[1]))
    ui.setup_ui(config_setter)
    config_setter.show()
    sys.exit(app.exec_())
