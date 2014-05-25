#!/usr/bin/kivy
# fileencoding=utf-8
__author__ = 'shibasay'

# written using copy & pasting and modifying code shown in below URL
# http://kivy.org/docs/api-kivy.uix.filechooser.html

from kivy.app import App
from kivy.uix.floatlayout import FloatLayout
from kivy.uix.gridlayout import GridLayout
from kivy.factory import Factory
from kivy.properties import ObjectProperty
from kivy.uix.popup import Popup
from kivy.uix.scatter import Scatter
from kivy.uix.image import Image
import os
import tempfile
import shutil

import bed_mod

class PreviewWindow(FloatLayout):
    def __init__(self, **kwargs):
        super(PreviewWindow, self).__init__()

        #setup Layouts
        client_frame      = GridLayout( size_hint = (1, 1), cols = 2, rows = 2)
        self.preview_xy_frame  = GridLayout( size_hint = (1, 1), cols = 1, rows = 1)
        self.preview_yz_frame  = GridLayout( size_hint = (1, 1), cols = 1, rows = 1)
        self.preview_zx_frame  = GridLayout( size_hint = (1, 1), cols = 1, rows = 1)
        self.add_widget(client_frame)
        client_frame.add_widget(self.preview_xy_frame)
        client_frame.add_widget(self.preview_yz_frame)
        client_frame.add_widget(self.preview_zx_frame)

    def showPNGonXY(self, path):
        scatter = Scatter()
        img = Image(source=path, allow_stretch=True)
        scatter.add_widget(img)
        self.preview_xy_frame.add_widget(scatter, 1)

    def showPNGonYZ(self, path):
        scatter = Scatter()
        img = Image(source=path, allow_stretch=True)
        scatter.add_widget(img)
        self.preview_yz_frame.add_widget(scatter, 1)
    def showPNGonZX(self, path):
        scatter = Scatter()
        img = Image(source=path, allow_stretch=True)
        scatter.add_widget(img)
        self.preview_zx_frame.add_widget(scatter, 1)


class LoadDialog(FloatLayout):
    load = ObjectProperty(None)
    cancel = ObjectProperty(None)


class SaveDialog(FloatLayout):
    save = ObjectProperty(None)
    text_input = ObjectProperty(None)
    cancel = ObjectProperty(None)


class Root(FloatLayout):
    loadfile = ObjectProperty(None)
    savefile = ObjectProperty(None)

    def __init__(self, **kwargs):
        super(Root, self).__init__()

        self.tmp_dir = tempfile.mkdtemp()

        #setup Layouts
        client_frame      = GridLayout( size_hint = (1, 1), cols = 2, rows = 1)
        self.preview_input     = PreviewWindow( size_hint = (1, 1))
        self.preview_output    = PreviewWindow( size_hint = (1, 1))
        self.add_widget(client_frame)
        client_frame.add_widget(self.preview_input)
        client_frame.add_widget(self.preview_output)

    def __del__(self):
        try:
            shutil.rmtree(self.tmp_dir)  # delete directory
        except OSError as exc:
            if exc.errno != errno.ENOENT:  # ENOENT - no such file or directory
                raise  # re-raise exception
        super(Root, self).__del__()

    def dismiss_popup(self):
        self._popup.dismiss()

    def show_load(self):
        content = LoadDialog(load=self.load, cancel=self.dismiss_popup)
        self._popup = Popup(title="Load file", content=content, size_hint=(0.9, 0.9))
        self._popup.open()

    def show_save(self):
        content = SaveDialog(save=self.save, cancel=self.dismiss_popup)
        self._popup = Popup(title="Save file", content=content, size_hint=(0.9, 0.9))
        self._popup.open()

    def load(self, path, filename):
        inpathname = os.path.join(path, filename[0])
        self.bed_input = bed_mod.bed_read(inpathname)
        self.bed_input_name = os.path.basename(inpathname)
        tmpdir = self.tmp_dir
        basename = self.bed_input_name
        self.preview_input.showPNGonXY(self.bed_input.getPreviewXY_path(basename, tmpdir))
        self.preview_input.showPNGonYZ(self.bed_input.getPreviewYZ_path(basename, tmpdir))
        self.preview_input.showPNGonZX(self.bed_input.getPreviewZX_path(basename, tmpdir))

        self.dismiss_popup()

    def save(self, path, filename):
        with open(os.path.join(path, filename), 'w') as stream:
            stream.write(self.bed_input.printAll())

        self.dismiss_popup()


class BedMod(App):
    pass

Factory.register('Root', cls=Root)
Factory.register('LoadDialog', cls=LoadDialog)
Factory.register('SaveDialog', cls=SaveDialog)

if __name__ == '__main__':
    BedMod().run()

