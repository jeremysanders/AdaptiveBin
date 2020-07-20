# Adaptive Binning

This page contains the source code for the algorithm described in [Sanders & Fabian (2001)](http://adsabs.harvard.edu/abs/2001MNRAS.325..178S). The code bins images with a varying bin size, depending on what the count-rate is. The implementation of the algorithm is called AdaptiveBin.
News

*    2011-09-14: Version 0.2.1 includes fixes for recent C++ compilers
*    2006: New Contour Binning algorithm released, which follows surface brightness contours. Please see https://github.com/jeremysanders/contbin
*    2001-12-13: Version 0.2.0 released.
*    2001-12-11: Program for extracting and fitting spectra in bins available. Please email me. I should release this when I get it into a easily usable form.
*    2001-02-11: Version 0.1.2 released.
*    2001-02-08: Paper accepted by MNRAS. Updated version below.
*    2001-01-10: New version 0.1.1 released. Adds support for making sure bins are contiguous.

## Building

A simple Makefile is provided. Some minor editing may be required to get it to compile.

## Documentation (AdaptiveBin)

 AdaptiveBin is the main Adbin program. It takes one or more images and adaptively bins them. If one image is supplied, then the pixels are binned by fractional error on the intensity. If two or more images are supplied, then the pixels are fractional binned by error on the combined colour.
Usage

```
$ AdaptiveBin --help
Usage: AdaptiveBin [OPTIONS] file bg=count...
Adaptively bins a set of images
Written by Jeremy Sanders, 2000, 2001.

  -o, --out=FILE           set out file (def adbin_out.fits)
  -e, --error=FILE         set error out file (def adbin_err.fits)
  -n, --binmap=FILE        set binmap out file (def adbin_binmap.fits)
  -p, --pixel=FILE         set pixel out file (DISCOURAGED)
  -t, --threshold=VAL      set threshold fraction (def 10%)
  -v, --value=STR          set output value (eg count(0), ratio(1,2))
  -s, --subpix=INT         set subpixel positioning divisior (def. 1)
  -c, --contig             only allow contiguous regions
      --verbose            display more information
      --help               display this help message
  -V, --version            display the program version

Report bugs to <jss@ast.cam.ac.uk>
```

AdaptiveBin takes a list of input FITS files and (optional) background counts. If one image is given, the image is binned based on the fractional error of the intensity of that image. If more than one images are supplied, then the adaptive binning is based on the fractional error of the combined colour of that image.

The switches `--out`, `--error` and `--binmap` set the output FITS images for the output binned image, error map and bin map, respectively. By default the file names are `adbin_out.fits`, `adbin_err.fits` and `adbin_binmap.fits`.

The maximum fractional error is set using the `--threshold=0.xx` option. By default it is 0.1.

The type of binned image produced is set using the `--value` switch. This option is only useful when colour binning is being done. If this option is set to count(x), then the pixels in the output-image will contain the average count in band x for their bins (x is numbered from 0). The error map will show the fractional error on that count (note the bins are always produced using the error on the combined-colour). If `--value` is set to `ratio(x,y)`, then the output-image will show the average colour x/y in the bin, and the error map will show the error on that colour.

The `--subpix=x` switch enables sub-pixel positioning (it should be called sub-bin positioning). The results of this option usually aren't great, but feel free to experiment. This option takes an integer value greater than 1. The value specifies where the top corner of a bin can be placed. Normally bins can only be placed on a regular grid, with a grid spacing of the size of the bin. This switch allows bins to be placed at intervals x times that frequency. Bins are drawn in ascending fractional error.

The `--contig` option makes sure that bins form contiguous regions. Using this option prevents 'stranded bins' which are binned together with a lower intensity region, due to them not having enough counts to have an error less than or equal to the threshold. If a bin consists of two isolated regions, then it is split into two different bins. A region is isolated if it does not have any neighbouring pixels (including sharing corners) with a different region. This option slows down the program, but there is probably some room for optimisation of the code.

### Notes

*    Version >= 0.1.2: AdaptiveBin can expand its options and arguments from a file, instead of the command line. Using an argument of `@filename` will substitute the text in the file in as options. The file can contain comments (preceeded by the # character); quote signs must be escaped using a backslash character.
*    Brackets are changed by shells. The --value option needs to be surrounded by quotes. For example, `--value="ratio(0,1)"`
*    WCS information is copied (in a simple form) from the first input image to the output images (including error and pixel maps).
*    All input images must be of the same size.

### Examples

```
$ AdaptiveBin --threshold=0.05 --out=adbin_intensity_005.fits \
    --err=adbin_intensity_error_005.fits --binmap=binmap_005.fits \
    infile.fits bg=0.874
```

Bin infile.fits using intensity-adaptive binning with a fractional error threshold of 0.05 (5%). The output filename is adbin_intensity_005.fits, with error and binmap filenames adbin_intensity_error_0.05.fits and binmap_005.fits. A background of 0.874 counts per pixel is used to analyse the input file.

```
$ AdaptiveBin --threshold=0.2 --value="ratio(1,2)" \
    image0.fits image1.fits bg=2.3 image2.fits bg=4.12
```

Do combined colour adaptive binning with a fractional error threshold of 0.2. The input images are image0.fits (no background), image1.fits (background 2.3 counts per pixel) and image2.fits (background 4.12 counts per pixel). The output, error and pixel filenames are the defaults. A colour image is outputted in the output and error images showing the ratio of the second and third input images (image1.fits and image2.fits). 

## ABPixelCopy documentation

ABPixelCopy is a utility which takes a pixel-map generated from AdaptiveBin and applies it to another image, creating an intensity map with the same bins.
Usage

```
Usage: ABPixelCopy [OPTIONS] --image=in.fits --out=out.fits --err=err.fits
Bins image using pixel file, creating output and error file.
Does background subtraction.
Written by Jeremy Sanders, 2000, 2001.

  -i, --image=FILE         set input image file (req)
  -o, --out=FILE           set out file (req)
  -e, --err=FILE           set out error map (req)
  -n, --binmap=FILE        set binmap file (def adbin_binmap.fits)
  -p, --pixel=FILE         set input pixel file (DISCOURAGED)
  -b, --background=VAL     set background counts/pixel
      --verbose            display more information
      --help               display this help message
  -V, --version            display the program version

Report bugs to <jss@ast.cam.ac.uk>
```

The syntax is a little weird, as everything is a switch, but some switches are required (marked req above). `--image` specifies the input image filename. `--out` and `--err` specify the output and error image filenames, respectively. `--binmap` specifies the input bin map. --background sets the X-ray background in counts per pixel to subtract from the input image and use to generate the error-map.

### Notes

*    Version >= 0.1.2: ABPixelCopy can expand its options and arguments from a file, instead of the command line. Using an argument of @filename will substitute the text in the file in as options. The file can contain comments (preceeded by the # character); quote signs must be escaped using a backslash character.
*    Copies WCS information from the input image to the output image.
*    Does not do colour binning.

### Examples

```
# ABPixelCopy --image=infile.fits --out=outfile.fits \
    --err=outfile_err.fits --binmap=binmap_050.fits -b 8.3
```

Take the bins in binmap_050.fits and apply them to infile.fits. The output image is outfile.fits, with the output error image outfile_err.fits. Assume a background of 8.3 counts per pixel for the input image.

