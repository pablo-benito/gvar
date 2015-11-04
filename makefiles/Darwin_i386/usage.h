const char * USAGE = { "VERBOSE     -h more verbose description of the following options\nINSTRUMENT  -i <instrument>  	       default: none\n                [none imager sounder iIcse \n		 sIcse iSpacelookStats iNav sNav]\n FILETYPE   -f <file>                  default: none\n                [none bin jpeg tiff tifftile tiffstrip]\n TIMESTAMP  -T <stamp>                 default: fromFilename\n                [fromFilename framestart region]\n START      -s <VISpixel,VISline>      default: from gvar\n END        -e <VISpixel,VISline>      default: from gvar\n REGION     -r <CentralVISPixel,CentralVISLine,Pixels,Lines>  \n				       default: from gvar\n REGION     -R <CentralVISLongitude,CentralVISLatitude,Pixels,Lines>  \n				       default: from gvar\n AREA       -A <[+/-] Mbytes>          default: inactive \n LINKNAME   -l toggle link to latest   default: none\n CHANNEL    -c <channel>               default: all\n                [all 1 2 .. ]\n  WORDTYPE  -w <wordtype>              default: none\n                [none char uchar short ushort int uint float double]\n  UNITS     -u <units>                 default: counts\n                [counts radiance     \n                 tstar tscene modeA    (for IR ) \n                 albedo albedo2]       (for VIS)\n  GAIN      -g <gain> {y=((x-b)/g)^G}  default: 1.0\n  BIAS      -b <bias>                  default: 0.0\n  GAMMA     -G <gamma>                 default: 1.0\n  XSCALE    -x <pixels>                default: 1.0\n  YSCALE    -y <lines>                 default: 1.0\n  DIRECTORY -p <output_path>           default: ./   \n  HOURS     -H <start,end>             default: 0,23\n                GMT hours\n  MAP       -M <map>                   default: none  \n                  [none grid1 grid2 grid12 ] \n  MAPVALUE  -m <counts>                default: 1023 \n  COMPRESS  -C <scheme|QUALITY>        default: none|75\n		for tiff:\n                [cittrle ccittfax3 ccittrlew lzw\n                 next packbits thunderscan pixarfilm]\n		for jpeg: QUALITY [0 - 100] \nSNOOZE      -S <directory>             default: off \nDEBUG       -D <mode>                  default: none \n               [none header lineDoc block0time SAD CRC parity]\nCOUNTER     -K toggle counter off      default: on\nCMDFILE     -F <fileName>              default: none\nXECUTE      -X <scriptName>            default: none\nPRIORITY    -P                         default: off\n\n"  }; 