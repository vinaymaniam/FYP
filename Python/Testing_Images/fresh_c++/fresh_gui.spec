# -*- mode: python -*-

# $~ make
# $~ pyinstaller fresh_gui.spec

block_cipher = None


a = Analysis(['fresh_gui.py'],
             pathex=['/Users/Matteo/Documents/MATLAB/upsampling/hires'],
             binaries=None,
             datas=[('fresh','.')],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='fresh_gui',
          debug=False,
          strip=False,
          upx=True,
          console=False )
app = BUNDLE(exe,
             name='fresh_gui.app',
             icon=None,
             bundle_identifier=None)
