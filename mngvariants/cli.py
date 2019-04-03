from . import arg_parser, variants

def main():
    args = arg_parser.parse()
    print('UUID: {}, type: {}'.format(args.uuid, type(args.uuid)))
    print('Samples: {}, type: {}'.format(args.samples, type(args.samples)))
    print('Reference: {}, type: {}'.format(args.reference, type(args.reference)))
    print('Workspace: {}, type: {}'.format(args.workspace, type(args.workspace)))

    variants.main(args)

if __name__ == '__main__':
    main()
