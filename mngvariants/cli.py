from . import arg_parser, variants


def main():
    args = arg_parser.parse()
    variants.main(args)


if __name__ == '__main__':
    main()
